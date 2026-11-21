"""Render source collision geometry as a visual check, without dynamics."""
import struct
import zlib
from pathlib import Path
import compare as c


def png(path, image):
    h, w, _ = image.shape
    def chunk(tag, data):
        return struct.pack('!I', len(data))+tag+data+struct.pack('!I', zlib.crc32(tag+data)&0xffffffff)
    raw = b''.join(b'\x00'+row.tobytes() for row in image)
    path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR', struct.pack('!2I5B',w,h,8,6,0,0,0))+
                     chunk(b'IDAT', zlib.compress(raw))+chunk(b'IEND', b''))


def main():
    c.RENDER = True
    c.pb.connect(c.pb.DIRECT)
    images = []
    for model, target, distance in [('franka_panda/panda.urdf',[0,0,.4],1.5), ('a1/a1.urdf',[0,0,.25],1.0)]:
        c.pb.resetSimulation()
        shapes, joints, roots, _, _ = c.load_robot(model)
        q = {j['name']:(j['lower']+j['upper'])/2 for j in joints if j['kind']!='fixed'}
        world = c.poses(joints, roots, q, c.transform([0,0,.3 if 'a1' in model else 0]))
        for shape in shapes:
            shape.update(world[shape.link] @ shape.local)
        view = c.pb.computeViewMatrixFromYawPitchRoll(target,distance,45,-25,0,2)
        proj = c.pb.computeProjectionMatrixFOV(50,1,.01,10)
        _,_,pixels,_,_ = c.pb.getCameraImage(600,600,viewMatrix=view,projectionMatrix=proj,renderer=c.pb.ER_TINY_RENDERER)
        images.append(c.np.array(pixels,dtype=c.np.uint8).reshape(600,600,4))
    png(Path(__file__).with_name('models.png'), c.np.concatenate(images,axis=1))
    c.pb.disconnect()


if __name__ == '__main__':
    main()
