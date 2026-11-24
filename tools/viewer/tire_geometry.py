"""Closed tire meshes with a smooth, compression-driven flat contact patch."""
import meshcat.geometry as geometry
import numpy as np


def smooth_clearance(clearance, blend):
    """Flatten penetration and join the undeformed surface with C2 continuity."""
    clearance = np.asarray(clearance, dtype=float)
    t = np.clip(clearance / blend, 0, 1)
    weight = t**3 * (10 + t * (-15 + 6 * t))
    return np.maximum(0, clearance) * weight


def tire_support(radius, rim_radius, width, pose, normal):
    """Lowest torus support point, matching the C++ tire contact geometry."""
    axle = pose[:3, 1]
    parallel = np.dot(normal, axle) * axle
    radial = normal - parallel
    major, minor, half_width = (radius + rim_radius) / 2, (radius - rim_radius) / 2, width / 2
    support = np.sqrt(minor**2 * np.dot(radial, radial) + half_width**2 * np.dot(parallel, parallel))
    arm = -(minor**2 * radial + half_width**2 * parallel) / support
    if np.linalg.norm(radial) > 1e-12:
        arm -= major * radial / np.linalg.norm(radial)
    return pose[:3, 3] + arm


def tire_mesh(radius, rim_radius, width, height=None, segments=96, *, pose=None,
              ground_normal=(0, 0, 1), ground_height=0, sections=16):
    """Closed elliptic torus deformed along the world ground normal.

    Returns world-space vertices when pose is supplied. The height-only form
    returns hub-relative vertices for compatibility with upright callers.
    """
    normal = np.asarray(ground_normal, dtype=float)
    if not (0 < rim_radius < radius and width > 0 and np.isfinite([radius, rim_radius, width, ground_height]).all()
            and np.isfinite(normal).all() and abs(np.linalg.norm(normal) - 1) < 1e-8):
        raise ValueError("Invalid tire geometry or ground plane")
    local_output = pose is None
    if local_output:
        pose = np.eye(4)
        pose[2, 3] = height
    theta = np.arange(segments) * 2 * np.pi / segments
    phi = np.arange(sections) * 2 * np.pi / sections
    major, minor = (radius + rim_radius) / 2, (radius - rim_radius) / 2
    radial = major + minor * np.cos(phi)
    vertices = np.empty((segments, sections, 3))
    vertices[:, :, 0] = np.cos(theta)[:, None] * radial
    vertices[:, :, 1] = width / 2 * np.sin(phi)[None, :]
    vertices[:, :, 2] = np.sin(theta)[:, None] * radial
    vertices = vertices.reshape(-1, 3) @ pose[:3, :3].T + pose[:3, 3]
    compression = max(0., ground_height - np.dot(normal, tire_support(radius, rim_radius, width, pose, normal)))
    if compression > 0:
        clearance = vertices @ normal - ground_height
        blend = min(.4 * minor, .5 * compression)
        deformed = smooth_clearance(clearance, blend)
        vertices += (deformed - clearance)[:, None] * normal
    faces = []
    for i in range(segments):
        for j in range(sections):
            a, b = i * sections + j, ((i + 1) % segments) * sections + j
            c = ((i + 1) % segments) * sections + (j + 1) % sections
            d = i * sections + (j + 1) % sections
            faces.extend(((a, c, b), (a, d, c)))
    if local_output:
        vertices -= pose[:3, 3]
    return vertices, np.asarray(faces)


def populate_tires(root, body_nodes, scene, animation, matrix):
    """Store sampled tire shapes with visibility tracks for portable Meshcat playback."""
    for tire_index, tire in enumerate(scene.get("tires", [])):
        body = tire["body"]
        radius, rim, width = tire["radius"], tire["rim_radius"], tire["width"]
        wheel = body_nodes[body]
        wheel["rim"].set_object(geometry.Cylinder(width * .3, rim), geometry.MeshLambertMaterial(color=0xAAB7C4))
        for spoke in range(3):
            node = wheel[f"spoke_{spoke}"]
            node.set_object(geometry.Box([2 * rim * .9, .012, .018]), geometry.MeshLambertMaterial(color=0x23364A))
            angle = spoke * np.pi / 3
            transform = np.eye(4)
            transform[:3, :3] = [[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]]
            transform[1, 3] = -width * .16
            node.set_transform(transform)
        material = geometry.MeshLambertMaterial(color=0x30343B)
        for index, frame in enumerate(scene["frames"]):
            pose = matrix(frame["poses"][body])
            vertices, faces = tire_mesh(radius, rim, width, pose=pose,
                                        ground_normal=tire.get("ground_normal", [0, 0, 1]),
                                        ground_height=tire["ground_height"])
            node = root[f"tire_{tire_index}"][f"sample_{index}"]
            node.set_object(geometry.TriangularMeshGeometry(vertices, faces), material)
            transform = np.eye(4)
            node.set_transform(transform)
            normal = np.asarray(tire.get("ground_normal", [0, 0, 1]), dtype=float)
            lowest = tire_support(radius, rim, width, pose, normal)
            gap = np.dot(normal, lowest) - tire["ground_height"]
            if gap < 0:
                patch = node["patch_center"]
                patch.set_object(geometry.Sphere(.009), geometry.MeshBasicMaterial(color=0xFF7733, transparent=True, depthTest=False, depthWrite=False))
                transform = np.eye(4)
                transform[:3, 3] = lowest - gap * normal
                patch.set_transform(transform)
            node.set_property("visible", index == 0)
            with animation.at_frame(node, 0) as animated:
                animated.set_property("visible", "bool", index == 0)
            if index:
                with animation.at_frame(node, index) as animated:
                    animated.set_property("visible", "bool", True)
            if index + 1 < len(scene["frames"]):
                with animation.at_frame(node, index + 1) as animated:
                    animated.set_property("visible", "bool", False)
