### Download data from Synapse
################################################################################

### Install and load necessary packages
# Install synapseclient python package
#reticulate::py_install("synapseclient", envname = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/env/basic/virtualenvs/r-reticulate")

### Load synapseclient and connect to synapse to download data
# See https://help.synapse.org/docs/Downloading-Data-Programmatically.2003796248.html
# Load synapseclient
synapseclient = reticulate::import("synapseclient")
# Connect to synapse (https://www.synapse.org/)
syn = synapseclient$Synapse()
# Log in
# you need to be registered and have generated a Personal Access Tokens with download rights
# See https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
# The token and default download destination / cache has to be set in a .synapseConfig file in the home directory
# See https://help.synapse.org/docs/Client-Configuration.1985446156.html
syn$login(authToken='eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIl0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTc0MDQ5MDYzMywiaWF0IjoxNzQwNDkwNjMzLCJqdGkiOiIxNjg0OSIsInN1YiI6IjM1MzIyNDkifQ.IF9NKtT-R4F1iX3Zma3_j1Zmb3_jG6jQ7fweucUEfnt7LZAV3ZsI753PwEDJFBk2lZzWpiwZVEZmbw8vIgZa2jYWyX3pnLnpIHmK3YTze4Pk-xtuqdKK7dtJLvdJZiWoxxRnXolMmEU0h9u-FYz-0oK7Q4BVfAHGghDsVUWiEDst_nqs0XSdP9IKdCU7BndnW0pASvocU5ZdaqHObxV9TSuymDjc5hv6XSn5dytMJxsJ9CcRDIlXPpnTmseGT62ANbL2LoPzBFmZS0eUHeqkZ13s9e76oU7bFCyI1SL8aGSQzYNyL_C1wyesZqzssN3XPkE_ATYECbZC2-putFBm-g')
# If token can be obtained form .synapseConfig
#syn$login()
# Specifiy download location and get data
destination_folder = "/mnt/ngsnfs/single_cell_dev/scRNAseq_processing/sc_analysis/references"
entity = syn$get(entity = "syn52943357", downloadLocation=destination_folder)

