from os import listdir
from os.path import isfile, join

import cv2
import numpy

mypath='build/images'
images = []
for n in range(0, 400):
    print(join(mypath, str(n) + ".png"))
    img = cv2.imread( join(mypath,str(n) + ".png") )
    blur = cv2.blur(img,(40,40))
    blur[blur < 123] = 0
    blur[blur >= 123] = 255
    images.append(blur)
fourcc = cv2.VideoWriter_fourcc('m','p','4','v')
writer = cv2.VideoWriter("water.mov", fourcc, 10, (images[0].shape[1], images[0].shape[0]))

for image in images:
    writer.write(image)
writer.release()
