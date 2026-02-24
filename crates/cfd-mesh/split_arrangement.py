import os

path = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-mesh\src\application\csg\arrangement.rs"
with open(path, "r", encoding="utf-8") as f:
    lines = f.readlines()

output_dir = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-mesh\src\application\csg\arrangement"
os.makedirs(output_dir, exist_ok=True)

header = "".join(lines[0:58])
classify_code = "".join(lines[58:241])
propagate_code = "".join(lines[243:589])
main_impl = "".join(lines[590:1650])

with open(os.path.join(output_dir, "classify.rs"), "w", encoding="utf-8") as f:
    f.write("use crate::domain::core::scalar::{Point3r, Real, Vector3r};\n")
    f.write("use crate::infrastructure::storage::face_store::FaceData;\n")
    f.write("use crate::infrastructure::storage::vertex_pool::VertexPool;\n")
    f.write(classify_code)

with open(os.path.join(output_dir, "propagate.rs"), "w", encoding="utf-8") as f:
    f.write("use crate::application::welding::snap::SnapSegment;\n")
    f.write("use crate::domain::core::scalar::{Point3r, Real, Vector3r};\n")
    f.write("use crate::infrastructure::storage::face_store::FaceData;\n")
    f.write("use crate::infrastructure::storage::vertex_pool::VertexPool;\n")
    f.write("use std::collections::HashMap;\n")
    f.write(propagate_code)

with open(os.path.join(output_dir, "mod.rs"), "w", encoding="utf-8") as f:
    f.write(header)
    f.write("pub mod classify;\n")
    f.write("pub mod propagate;\n")
    f.write("#[cfg(test)]\n")
    f.write("pub mod tests;\n\n")
    f.write("use classify::*;\n")
    f.write("use propagate::*;\n\n")
    f.write(main_impl)

with open(os.path.join(output_dir, "tests.rs"), "w", encoding="utf-8") as f:
    f.write("use super::*;\n")
    unwrapped_tests = lines[1653:-1]
    f.writelines(unwrapped_tests)

os.remove(path)
