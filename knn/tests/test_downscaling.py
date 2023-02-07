import os
import shutil
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from osgeo import osr
import numpy
import pandas
import pygeoprocessing
from shapely.geometry import Polygon
import xarray


class TestKNN(unittest.TestCase):
    """Tests knn.py."""

    def setUp(self):
        self.workspace_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workspace_dir)

    def test_shift_longitude_error_on_missing_dimension(self):
        """Test shift longitude from 0-360."""
        from .. import knn

        precip = 10 * numpy.random.rand(2, 2)
        lon = [[-99.83, -99.32], [-99.79, -99.23]]
        lat = [[42.25, 42.21], [42.63, 42.59]]
        dataset = xarray.Dataset(
            {
                'pr': (['x', 'y'], precip),
            },
            coords={
                'long': (['x', 'y'], lon),  # deliberately mislabeled coord.
                'lat': (['x', 'y'], lat),
            }
        )

        with self.assertRaises(ValueError) as cm:
            _ = knn.shift_longitude_from_360(dataset)
        self.assertTrue(
            'Expected dimension "lon" but found coordinates' in str(cm.exception))

    def test_shift_longitude_correct_given_180(self):
        """Test shift longitude is correct given coords from -180-180."""
        from .. import knn

        precip = 10 * numpy.random.rand(2, 2)
        lon = [[-180.0, -5.0], [0.0, 179.9]]
        lat = [[42.25, 42.21], [42.63, 42.59]]
        dataset = xarray.Dataset(
            {
                'pr': (['x', 'y'], precip),
            },
            coords={
                'lon': (['x', 'y'], lon),
                'lat': (['x', 'y'], lat),
            }
        )

        shifted_dataset = knn.shift_longitude_from_360(dataset)
        numpy.testing.assert_array_almost_equal(
            shifted_dataset.lon.to_numpy(), lon)

    def test_tri_state_joint_probability(self):
        """Test tri_state_joint_probability."""
        from .. import knn

        timeseries = numpy.array(
            [0, 0, 2, 5, 2, 1, 0, 0, 9, 8, 0])
        lower_bound = 1.5
        upper_bound = 7.5

        expected_matrix = numpy.array([
            [0.3, 0.1, 0.1],
            [0.1, 0.2, 0.0],
            [0.1, 0.0, 0.1]
        ])

        actual_matrix = knn.tri_state_joint_probability(
            timeseries, lower_bound, upper_bound)

        numpy.testing.assert_array_equal(actual_matrix, expected_matrix)

    def test_slice_dates_around_dayofyear(self):
        """Test slicing dates from an xarray.Dataset time coordinate"""
        from .. import knn

        dates = ('1979-01-01', '1981-12-31')
        dates_index = knn.date_range_no_leap(*dates)
        dataset = xarray.Dataset({'time': dates_index})
        month = 10
        day = 16
        near_window = 4
        expected_dates = [
            '1979-10-12', '1979-10-13', '1979-10-14', '1979-10-15',
            '1979-10-16', '1979-10-17', '1979-10-18', '1979-10-19', '1979-10-20',
            '1980-10-12', '1980-10-13', '1980-10-14', '1980-10-15',
            '1980-10-16', '1980-10-17', '1980-10-18', '1980-10-19', '1980-10-20',
            '1981-10-12', '1981-10-13', '1981-10-14', '1981-10-15',
            '1981-10-16', '1981-10-17', '1981-10-18', '1981-10-19', '1981-10-20',
        ]
        idx = knn.slice_dates_around_dayofyear(
            dataset.time, month, day, near_window)
        self.assertListEqual(
            [d.strftime(format='%Y-%m-%d') for d in dates_index[idx]],
            expected_dates)

    def test_state_transition_series(self):
        """Test create series of state transitions."""
        from .. import knn

        array = numpy.array([0, 0, 0, 5, 2, 1, 0, 0, 9, 8, 0])
        expected_array = [
            'AA', 'AA', 'AB', 'BA', 'AA',
            'AA', 'AA', 'AC', 'CC', 'CA']
        lower_bound = 3
        upper_bound = 7
        transitions = knn.state_transition_series(
            array, lower_bound, upper_bound)
        numpy.testing.assert_array_equal(transitions, expected_array)

    def test_jp_matrix_from_transitions_sequence(self):
        """Test joint probablity matrix from sequence of transitions."""
        from .. import knn

        n = 900
        reference_dates = ('1980-01-01', '1989-12-31')
        ref_dates = knn.date_range_no_leap(*reference_dates)
        ref_dates = ref_dates[:n]
        dataset = xarray.Dataset({'time': ref_dates})
        month = 6
        day = 20
        near_window = 4

        sample_transitions = [
            'AA', 'AB', 'AC', 'BA', 'BB', 'BC', 'CA', 'CB', 'AA']
        transitions = numpy.empty(ref_dates.shape, dtype='U2')
        transitions[:] = 100 * sample_transitions

        # Our sample_transitions pattern repeats every 9 days, and a full
        # window sequence is 9 days, so the expected frequencies of
        # transitions in the computed matrix match the frequencies of
        # occurrence in sample_transitions:
        expected_matrix = numpy.array([
            [0.222222, 0.111111, 0.111111],
            [0.111111, 0.111111, 0.111111],
            [0.111111, 0.111111, 0.0]
        ])
        jp_matrix = knn.jp_matrix_from_transitions_sequence(
            dataset.time, transitions, month, day, near_window)
        numpy.testing.assert_array_almost_equal(jp_matrix, expected_matrix)

    @patch('knn.knn.GCS_PROTOCOL', 'file:///')
    def test_execute(self):
        """Test the whole execute method."""

        from knn import knn
        thing = knn.GCSFS

        observed_dataset_path = os.path.join(self.workspace_dir, 'obs.nc')
        historical_gcm_path = os.path.join(
            self.workspace_dir, 'historical_gcm.zarr')
        future_gcm_path = os.path.join(
            self.workspace_dir, 'future_gcm.zarr')
        thing.glob = MagicMock(
            side_effect=[[historical_gcm_path], [future_gcm_path]])

        aoi_path = os.path.join(self.workspace_dir, 'aoi.geojson')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        wkt = srs.ExportToWkt()
        aoi_geometries = [Polygon([
            (1, 1), (1, 5), (5, 5), (5, 1), (1, 1)])]
        pygeoprocessing.shapely_geometry_to_vector(
            aoi_geometries, aoi_path, wkt, 'GeoJSON')

        start = '1980-01-01'
        end = '1989-12-31'
        dates = pandas.date_range(start, end, freq='D')
        lons = [0.0, 2.0, 4.0, 6.0]
        lats = [0.0, 2.0, 4.0, 6.0]

        precip_array = numpy.random.negative_binomial(
            n=1, p=0.4, size=(len(dates), len(lons), len(lats)))
        da = xarray.DataArray(
            precip_array,
            coords={
                'time': dates,
                'lon': lons,
                'lat': lats
            },
            dims=('time', 'lon', 'lat')
        )
        precip_dataset = xarray.Dataset({knn.MSWEP_VAR: da})
        precip_dataset.to_netcdf(observed_dataset_path)

        historic_gcm_array = numpy.random.negative_binomial(
            n=1, p=0.4, size=(len(dates), len(lons), len(lats)))
        historic_da = xarray.DataArray(
            historic_gcm_array,
            coords={
                'time': dates,
                'lon': lons,
                'lat': lats
            },
            dims=('time', 'lon', 'lat')
        )
        historic_gcm_dataset = xarray.Dataset({knn.GCM_PRECIP_VAR: historic_da})
        historic_gcm_dataset.to_zarr(historical_gcm_path)

        forecast_start = '2000-01-01'
        forecast_end = '2001-12-31'
        forecast_dates = pandas.date_range(
            forecast_start, forecast_end, freq='D')
        precip_array = numpy.random.negative_binomial(
            n=3, p=0.4, size=(len(forecast_dates), len(lons), len(lats)))
        forecast_da = xarray.DataArray(
            precip_array,
            coords={
                'time': forecast_dates,
                'lon': lons,
                'lat': lats
            },
            dims=('time', 'lon', 'lat')
        )
        forecast_dataset = xarray.Dataset({knn.GCM_PRECIP_VAR: forecast_da})
        forecast_dataset.to_zarr(future_gcm_path)

        args = {
            'reference_period_dates': (start, end),
            'prediction_dates': (forecast_start, forecast_end),
            'hindcast': True,
            'gcm_experiment_list': [knn.GCM_EXPERIMENT_LIST[0]],
            'gcm_model_list': [knn.MODEL_LIST[0]],
            'upper_precip_percentile': 75,
            'lower_precip_threshold': 1,  # millimeter
            'aoi_path': aoi_path,
            'observed_dataset_path': observed_dataset_path,
            'workspace_dir': self.workspace_dir,
        }

        knn.execute(args)

