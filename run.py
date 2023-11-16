#!/usr/bin/env python3

import subprocess
import argparse
from PIL import Image

# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("--method", type=str)
ap.add_argument("--lambd", type=float)
ap.add_argument("--denoise", type=int)
ap.add_argument("--sigma", type=int)
args = ap.parse_args()

add_noise = 1 if not args.denoise else 0

methods = ['111', '211', '221', '-111',
               '-1-11', '-121', '2-11', '110', '-110']
pqrSet = [(1, 1, 1), (2, 1, 1), (2, 2, 1), (-1, 1, 1),
            (-1, -1, 1), (-1, 2, 1), (2, -1, 1), (1, 1, 0), (-1, 1, 0)]
ids = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
methodsids = [{'method' : m, 'id' : mid}
            for m, mid in zip(methods, ids)]

files = ['input_0']

#resize input to maximum height=400pixels and maximum width=400pixels
m = 'input_0'
img = Image.open(f'{m}.png')
(dx, dy) = img.size
if dy > 400:
    img.resize((int(dx*400/dy), 400))
    (dx, dy) = img.size
    if dx > 400:
        img = img.resize((400, int(dy*400/dx)))
    img.save(f'{m}.png')

outname = 'denoised'
diffname = 'diff'

if args.denoise:
    #denoise / cartoon+texture
    #p=[]
    for m in methodsids:
        if str(m['method']) == args.method:
            i = int(m['id'])
           
            
            files += [outname]
            p = ['denoisingPDHG_ipol', 'input_0.png', outname + '.png', str(args.lambd), str(pqrSet[i][0]), str(pqrSet[i][1]), str(pqrSet[i][2])]
            subprocess.run(p)

    #compute difference denoised-original
    sigma = 5
    p = []
    for m in methodsids:
        if str(m['method']) == args.method:
            i = int(m['id'])

            #self.cfg['diff']['%i'%i] = diffname
            files += [diffname]

            p = ['imdiff_ipol', 'input_0.png', outname + '.png', diffname + '.png', str(sigma)]
            subprocess.run(p)

if add_noise:
    #add noise
    files += ['noisy']
    p = ['addnoise_ipol', 'input_0.png', 'noisy.png', str(args.sigma)]
    subprocess.run(p)

    #denoise
    #p=[]

    for m in methodsids:
        #if str(self.cfg['param'][m['method']]) == 'True':
        if str(m['method']) == args.method:
            i = int(m['id'])

            files += [outname]

            p = ['denoisingPDHG_ipol', 'noisy.png', outname + '.png', str(args.lambd), 
                                str(pqrSet[i][0]), str(pqrSet[i][1]), str(pqrSet[i][2])]
            subprocess.run(p)

    #compute difference denoised-original
    p = []
    for m in methodsids:
        #if str(self.cfg['param'][m['method']]) == 'True':
        if str(m['method']) == args.method:
            i = int(m['id'])

            files += [diffname]

            with open(f'{diffname}_rmse.txt', 'w') as file:
                p = subprocess.run(['imdiff_ipol', 'input_0.png', outname + '.png', diffname + '.png', str(args.sigma)], stdout=file)

    for m in methodsids:
        #if str(self.cfg['param'][m['method']]) == 'True':
        if str(m['method']) == args.method:
            i = int(m['id'])


# Resize for visualization (always zoom by at least 2x)
(sizeX, sizeY) = Image.open('input_0.png').size
zoomfactor = 1
if max(sizeX, sizeY) < 250:
    zoomfactor = 2
#zoomfactor = max(1, int(math.ceil(400.0/max(sizeX, sizeY))))

if zoomfactor > 1:
    (sizeX, sizeY) = (zoomfactor*sizeX, zoomfactor*sizeY)

    for filename in files:
        im = Image.open(filename + '.png')
        im.resize((sizeX, sizeY), method='nearest')
        im.save(filename + '_zoom.png')

