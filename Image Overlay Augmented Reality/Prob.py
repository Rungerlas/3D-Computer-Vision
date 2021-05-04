# -*- coding: utf-8 -*-
"""
2020年5月8日18:43:29

@author: Tu Kechao ，Ruohuan Xu
"""

from cv2 import cv2
import numpy as np
import math
import time
import matplotlib.pyplot as plt
import random
import os
import copy
import sift

from PCV.geometry import homography, camera
from pyexiv2 import Image
from pylab import *
from homework3 import CornersDetector

def perf_time(func):
    def wrap(*args):
        start = time.time()
        result = func(*args)
        cost = time.time() - start
        print("{} used {} s".format(func.__name__, cost))
        return result
    return wrap


def loadPics(picPath):
    img = cv2.imread(picPath,0)
    return img

def getFeaturePointsCoordinates(matrix):
    h,w = matrix.shape
    l = []
    for i in range(h):
        for j in range(w):
            if matrix[i,j] == True:
                l.append((j,i))
    # re = np.zeros((2,len(l)))
    # for k in range(len(l)):
    #     re[0,k] = l[k][0]
    #     re[1,k] = l[k][1]
    return l

def get_widvalue(I,c_coords1,c_coords2,wid=5):#获取角点周围点信息存进desc
    desc1 = []
    for i in range(len(c_coords1)):
        patch=I[c_coords1[i]-wid:c_coords1[i]+wid+1,c_coords2[i]-wid:c_coords2[i]+wid+1]
        desc1.append(patch)
    return desc1

def getPatchValue(Img,ListOfCoords,r=3):
    re = []
    for i in range(len(ListOfCoords)):
        x = int(ListOfCoords[i][1])
        y = int(ListOfCoords[i][0])
        patch = Img[x-r:x+r+1,y-r:y+r+1]
        re.append(patch)
    return re

def SADloop(Left,Right,ListOfCoords1,ListOfCoords2):
    Patch1 = getPatchValue(Left,ListOfCoords1)
    Patch2 = getPatchValue(Right,ListOfCoords2)
    matchList = []
    patchNum = len(Patch1)
    flag = 1
    for i in range(patchNum):
        if  i == patchNum * flag/10:
            print("完成了%f%%" % (flag*10))
            flag += 1
        LKernal = Patch1[i]
        sadTemp = []
        for j in range(len(Patch2)):
            RKernal = Patch2[j]
            sad = np.sum(np.abs(RKernal.astype(
                        np.double) - LKernal.astype(np.double)))
            sadTemp.append(sad)
        matchList.append(np.argmin(sadTemp))
    matchList = np.asarray(matchList)
    right = []
    left = ListOfCoords1
    for k in range(patchNum):
        # re[0,k] = ListOfCoords2[0,matchList[k]]
        # re[1,k] = ListOfCoords2[1,matchList[k]]
        right.append((int(ListOfCoords2[matchList[k]][0]),int(ListOfCoords2[matchList[k]][1])))
    print(len(left))
    print(len(right))
    return ListOfCoords1,right


def DrawImageCorresponding(Img1,Img2,Left,Right):
    h1,w1,_ = Img1.shape
    imgCombined = np.zeros((h1,2*w1,3),dtype=np.uint8)
    for i in range(h1):
        for j in range(w1):
            #draw left
            imgCombined[i,j] = Img1[i,j]
            #draw right
            imgCombined[i,j+w1] = Img2[i,j]
    for k in range(len(Left)):
        l = (Left[k][0],Left[k][1])
        r = (Right[k][0]+w1,Right[k][1])
        cv2.line(imgCombined,l,r,getRandomColor(),thickness=1)
    return imgCombined


def getRandomColor():
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255))

def cameraPoseFromHomography(H):
    H1 = H[:, 0]
    H2 = H[:, 1]
    H3 = np.cross(H1, H2)
    norm1 = np.linalg.norm(H1)
    norm2 = np.linalg.norm(H2)
    tnorm = (norm1 + norm2) / 2.0
    T = H[:, 2] / tnorm
    return np.mat([H1, H2, H3, T])


def getInternalCalibrationMatrix(path):
    #Using the given info
    # CMOS长宽信息
    w_c = 17.3
    h_c = 13.0
    # 读取图片的长宽信息
    img = cv2.imread(path)
    h = img.shape[0]
    w = img.shape[1]
    print(w, h)

    # 从EXIF中读取焦距
    i = Image(path)
    a, b = i.read_exif().get('Exif.Photo.FocalLength').split('/')
    fm = int(a) / int(b)
    f = w * fm / w_c

    # 计算内参信息
    K = np.zeros((3, 3))
    K[0][0] = f
    K[1][1] = f
    K[0][2] = w/2
    K[1][2] = h/2
    K[2][2] = 1
    return K


# def build3D(img1):
#     origin = [(1495,720),(1296,1063),(1738,1204),(1911,842)]
#     drawCube(img1,origin,500,500,(0,255,0))
    
#     plt.imshow(img1)
#     plt.show()

# def drawCube(img,origin,offsetX,offsetY,color,t=2,h=500):
#     p0 = (origin[0][0]+offsetX,origin[0][1]+offsetY)
#     p1 = (origin[1][0]+offsetX,origin[1][1]+offsetY)
#     p2 = (origin[2][0]+offsetX,origin[2][1]+offsetY)
#     p3 = (origin[3][0]+offsetX,origin[3][1]+offsetY)
#     cv2.line(img,p0,p1,color,thickness=t)
#     cv2.line(img,p1,p2,color,thickness=t)
#     cv2.line(img,p2,p3,color,thickness=t)
#     cv2.line(img,p3,p0,color,thickness=t)
#     # p0h = (origin[0][0]+offsetX,origin[0][1]+offsetY+h)
#     # p1h = (origin[1][0]+offsetX,origin[1][1]+offsetY+h)
#     # p2h = (origin[2][0]+offsetX,origin[2][1]+offsetY+h)
#     # p3h = (origin[3][0]+offsetX,origin[3][1]+offsetY+h)
#     # cv2.line(img,p0h,p1h,color,thickness=t)
#     # cv2.line(img,p1h,p2h,color,thickness=t)
#     # cv2.line(img,p2h,p3h,color,thickness=t)
#     # cv2.line(img,p3h,p0h,color,thickness=t)

def cube_points(c, wid):
    """ Creates a list of points for plotting
        a cube with plot. (the first 5 points are
        the bottom square, some sides repeated). """
    p = []
    # bottom
    p.append([c[0]-wid, c[1]-wid, c[2]-wid])
    p.append([c[0]-wid, c[1]+wid, c[2]-wid])
    p.append([c[0]+wid, c[1]+wid, c[2]-wid])
    p.append([c[0]+wid, c[1]-wid, c[2]-wid])
    p.append([c[0]-wid, c[1]-wid, c[2]-wid]) #same as first to close plot
    # top
    p.append([c[0]-wid, c[1]-wid, c[2]+wid])
    p.append([c[0]-wid, c[1]+wid, c[2]+wid])
    p.append([c[0]+wid, c[1]+wid, c[2]+wid])
    p.append([c[0]+wid, c[1]-wid, c[2]+wid])
    p.append([c[0]-wid, c[1]-wid, c[2]+wid]) #same as first to close plot
    # vertical sides
    p.append([c[0]-wid, c[1]-wid, c[2]+wid])
    p.append([c[0]-wid, c[1]+wid, c[2]+wid])
    p.append([c[0]-wid, c[1]+wid, c[2]-wid])
    p.append([c[0]+wid, c[1]+wid, c[2]-wid])
    p.append([c[0]+wid, c[1]+wid, c[2]+wid])
    p.append([c[0]+wid, c[1]-wid, c[2]+wid])
    p.append([c[0]+wid, c[1]-wid, c[2]-wid])
    return array(p).T

def projection_matrix(camera_parameters, homography):
    """
    From the camera calibration matrix and the estimated homography
    compute the 3D projection matrix
    """
    # Compute rotation along the x and y axis as well as the translation
    homography = homography * (-1)
    rot_and_transl = np.dot(np.linalg.inv(camera_parameters), homography)
    col_1 = rot_and_transl[:, 0]
    col_2 = rot_and_transl[:, 1]
    col_3 = rot_and_transl[:, 2]
    # normalise vectors
    l = math.sqrt(np.linalg.norm(col_1, 2) * np.linalg.norm(col_2, 2))
    rot_1 = col_1 / l
    rot_2 = col_2 / l
    translation = col_3 / l
    # compute the orthonormal basis
    c = rot_1 + rot_2
    p = np.cross(rot_1, rot_2)
    d = np.cross(c, p)
    rot_1 = np.dot(c / np.linalg.norm(c, 2) + d / np.linalg.norm(d, 2), 1 / math.sqrt(2))
    rot_2 = np.dot(c / np.linalg.norm(c, 2) - d / np.linalg.norm(d, 2), 1 / math.sqrt(2))
    rot_3 = np.cross(rot_1, rot_2)
    # finally, compute the 3D projection matrix from the model to the current frame
    projection = np.stack((rot_1, rot_2, rot_3, translation)).T
    return np.dot(camera_parameters, projection)

def render(img, obj, projection, model):
    vertices = obj.vertices
    scale_matrix = np.eye(3) * 3
    h, w = model.shape

    for face in obj.faces:
        face_vertices = face[0]
        points = np.array([vertices[vertex - 1] for vertex in face_vertices])
        points = np.dot(points, scale_matrix)
        # render model in the middle of the reference surface. To do so,
        # model points must be displaced
        points = np.array([[p[0] + w / 2, p[1] + h / 2, p[2]] for p in points])
        dst = cv2.perspectiveTransform(points.reshape(-1, 1, 3), projection)
        imgpts = np.int32(dst)
        cv2.fillConvexPoly(img, imgpts, (137, 27, 211))
    return img


def main(path1,path2,name,thre = 2000000):


    K = getInternalCalibrationMatrix(path1)

    img1 = cv2.imread(path1)
    img2 = cv2.imread(path2)
    img1H,corners1 = CornersDetector(path1,thre)
    img2H,corners2 = CornersDetector(path2,thre)
    cornersXY1 = getFeaturePointsCoordinates(corners1)
    cornersXY2 = getFeaturePointsCoordinates(corners2)
    points1,points2 = SADloop(img1,img2,cornersXY1,cornersXY2)
    # im = DrawImageCorresponding(img1,img2,points1,points2)

    src1 = np.float32(points1)
    src2 = np.float32(points2)
    H = cv2.findHomography(src2,src1,cv2.RANSAC,5.0)
    num, Rs, Ts, Ns  = cv2.decomposeHomographyMat(H[0], K)

    for i in range(num):
        print("R%d is " % i)
        print(Rs[i])
        print("T%d is "% i)
        print(Ts[i])

    model = cv2.imread("model.png", 0)
    print(K)
    P = projection_matrix(K,H[0])
    obj = OBJ(os.path.join(".\\fox.obj"), swapyz=False)
    im = render(img1, obj, P, model)

    # plt.subplot(2,1,1)
    # plt.imshow(img1H)
    # plt.subplot(2,1,2)
    # plt.imshow(img2H)
    # plt.imshow(im)
    # plt.show()
    cv2.imwrite(name,im)

class OBJ:
    def __init__(self, filename, swapyz=True):
        """Loads a Wavefront OBJ file. """
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []
        material = None
        for line in open(filename, "r"):
            if line.startswith('#'): continue
            values = line.split()
            if not values: continue
            if values[0] == 'v':
                v = list(map(float, values[1:4]))
                if swapyz:
                    v = v[0], v[2], v[1]
                self.vertices.append(v)
            elif values[0] == 'vn':
                v = list(map(float, values[1:4]))
                if swapyz:
                    v = v[0], v[2], v[1]
                self.normals.append(v)
            elif values[0] == 'vt':
                self.texcoords.append(map(float, values[1:3]))
            #elif values[0] in ('usemtl', 'usemat'):
                #material = values[1]
            #elif values[0] == 'mtllib':
                #self.mtl = MTL(values[1])
            elif values[0] == 'f':
                face = []
                texcoords = []
                norms = []
                for v in values[1:]:
                    w = v.split('/')
                    face.append(int(w[0]))
                    if len(w) >= 2 and len(w[1]) > 0:
                        texcoords.append(int(w[1]))
                    else:
                        texcoords.append(0)
                    if len(w) >= 3 and len(w[2]) > 0:
                        norms.append(int(w[2]))
                    else:
                        norms.append(0)
                #self.faces.append((face, norms, texcoords, material))
                self.faces.append((face, norms, texcoords))


if __name__ == "__main__":
    for i in range(13):
        n = 190+i
        fn1 = "Images\\P1070"+ str(n)+".JPG"
        fn2 = "Images\\P1070"+ str(n+1)+".JPG"
        name = "result"+str(n)+".jpg"
        main(fn1,fn2,name)