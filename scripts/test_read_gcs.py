# https://cloud.google.com/docs/authentication/provide-credentials-adc
# gcloud projects list
# gcloud config set project natcap-server
# gcloud auth application-default login

import pandas

df = pandas.read_csv(
    'gcs://natcap-climate-data/junk.txt',
    storage_options={
        'token': None, 'project': 'natcap-servers'
    })
