"""Orthographic camera and depth-buffer CAD rasterizer."""
import math
import numpy as np
from numba import njit

def normalized(v):
    p = np.asarray(v, dtype=float)
    return p / np.linalg.norm(p)


class Camera:
    def __init__(self, width, height, center, direction, horizontal, up=(0,-1,0)):
        self.width, self.height = width, height
        self.center = np.asarray(center)
        forward = normalized(direction)
        right = normalized(np.cross(forward, up))
        vertical = np.cross(right, forward)
        self.axes = np.vstack([right, -vertical, forward])
        self.scale = width/horizontal

    def project(self, xyz):
        projected = (np.asarray(xyz)-self.center) @ self.axes.T
        return projected * self.scale + np.array([self.width*.5,self.height*.5,0.0])


@njit
def rasterize(positions, ids, palette, w, h):
    zbuf = np.full((h,w), -1.e20, dtype=np.float32)
    labels = np.zeros((h,w), dtype=np.int32)
    picture = np.empty((h,w,3), dtype=np.uint8)
    picture[:,:,0], picture[:,:,1], picture[:,:,2] = 241,246,249
    for i in range(len(positions)):
        a,b,c = positions[i]
        x0 = max(0, int(math.floor(min(a[0],b[0],c[0]))))
        x1 = min(w-1, int(math.ceil(max(a[0],b[0],c[0]))))
        y0 = max(0, int(math.floor(min(a[1],b[1],c[1]))))
        y1 = min(h-1, int(math.ceil(max(a[1],b[1],c[1]))))
        if x1 < x0 or y1 < y0:
            continue
        u,v = b-a,c-a
        area = u[0]*v[1]-u[1]*v[0]
        if abs(area) < 0.00001:
            continue
        nx,ny,nz = u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],area
        if nz < 0:
            nx,ny,nz = -nx,-ny,-nz
        shade = max(0., (-.32*nx-.55*ny+.77*nz)/math.sqrt(nx*nx+ny*ny+nz*nz))
        light = .60 + .44*shade
        for y in range(y0,y1+1):
            dy = y+.5-a[1]
            for x in range(x0,x1+1):
                dx = x+.5-a[0]
                s = (dx*v[1]-dy*v[0])/area
                t = (dy*u[0]-dx*u[1])/area
                if s >= -1.e-6 and t >= -1.e-6 and s+t <= 1.000001:
                    z = a[2]+s*u[2]+t*v[2]
                    if z > zbuf[y,x]:
                        zbuf[y,x] = z
                        labels[y,x] = ids[i]
                        for ch in range(3):
                            picture[y,x,ch] = min(255., palette[ids[i],ch]*light)
    for y in range(h-1):
        for x in range(w-1):
            if labels[y,x] and (labels[y,x]!=labels[y,x+1] or labels[y,x]!=labels[y+1,x]):
                for ch in range(3):
                    picture[y,x,ch] = picture[y,x,ch]*.78
    return picture
