from ipykernel.comm import Comm
import requests

# Use comm to send a message from the kernel
xena_comm = Comm(target_name='xena_comm') #what is data?
session = requests.Session()
adapter = requests.adapters.HTTPAdapter(
    pool_connections=100,
    pool_maxsize=100)
session.mount('http://', adapter)
# This runs asynchronously, so errors aren't displayed in the notebook. Where are they logged?
def receiver(msg):
    query = msg['content']['data']['msg']
    resp = 'unset'
    if query['method'] == 'GET':
        resp = session.get(query['url'])
    if query['method'] == 'POST':
        resp = session.post(query['url'], data=query['body'], headers=query['headers'])
    xena_comm.send({'id': msg['content']['data']['id'], 'status': resp.status_code, 'response': resp.text})

# Add a callback for received messages.
@xena_comm.on_msg
def _recv(msg):
    receiver(msg)
