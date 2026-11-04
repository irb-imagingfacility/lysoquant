from ptparse import load
d = load('model.prototxt')
print(d['name'], len(d['layer']))
from collections import Counter
rows=[]
for L in d['layer']:
    inc = L.get('include')
    phase = inc[0]['phase'][0] if inc else '-'
    rows.append((L['name'][0], L['type'][0], phase))
print(Counter((t,p) for _,t,p in rows))
for L in d['layer']:
    inc = L.get('include'); phase = inc[0]['phase'][0] if inc else '-'
    if phase=='TRAIN': continue
    t=L['type'][0]
    extra=''
    for k in L:
        if k.endswith('_param') and k!='hdf5_data_param': extra+=f" {k}={ {a:b for a,b in L[k][0].items() if a!='weight_filler' and a!='bias_filler'} }"
    print(L['name'][0], t, L.get('bottom'), L.get('top'), extra)
