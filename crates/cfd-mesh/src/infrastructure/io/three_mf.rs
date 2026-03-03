//! 3MF (3D Manufacturing Format) export.
//!
//! 3MF is the modern replacement for STL in additive manufacturing workflows.
//! It is a ZIP archive containing an XML model file at `3D/3dmodel.model`.

use std::collections::HashMap;
use std::io::{Seek, Write};

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::mesh::IndexedMesh;

/// Write an [`IndexedMesh`] as a 3MF archive.
///
/// The output is a ZIP archive containing the required content types and
/// relationships files, plus the `3D/3dmodel.model` XML payload.
pub fn write_3mf<W: Write + Seek>(writer: W, mesh: &IndexedMesh) -> MeshResult<()> {
    let mut zip = zip::ZipWriter::new(writer);
    let options = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);

    // [Content_Types].xml
    zip.start_file("[Content_Types].xml", options)
        .map_err(io_err)?;
    zip.write_all(CONTENT_TYPES_XML.as_bytes())
        .map_err(MeshError::Io)?;

    // _rels/.rels
    zip.start_file("_rels/.rels", options).map_err(io_err)?;
    zip.write_all(RELS_XML.as_bytes()).map_err(MeshError::Io)?;

    // 3D/3dmodel.model
    zip.start_file("3D/3dmodel.model", options)
        .map_err(io_err)?;

    let model_xml = build_model_xml(mesh);
    zip.write_all(model_xml.as_bytes()).map_err(MeshError::Io)?;

    zip.finish().map_err(io_err)?;
    Ok(())
}

fn io_err(e: zip::result::ZipError) -> MeshError {
    MeshError::Other(format!("ZIP error: {e}"))
}

fn build_model_xml(mesh: &IndexedMesh) -> String {
    let mut xml = String::new();
    xml.push_str(MODEL_HEADER);

    // Vertices
    xml.push_str("      <vertices>\n");

    let mut id_to_idx: HashMap<crate::domain::core::index::VertexId, usize> = HashMap::new();
    let mut idx = 0usize;

    for (vid, vdata) in mesh.vertices.iter() {
        id_to_idx.insert(vid, idx);
        idx += 1;
        let p = &vdata.position;
        xml.push_str(&format!(
            "        <vertex x=\"{}\" y=\"{}\" z=\"{}\" />\n",
            p.x, p.y, p.z
        ));
    }
    xml.push_str("      </vertices>\n");

    // Triangles
    xml.push_str("      <triangles>\n");
    for (_fid, face) in mesh.faces.iter_enumerated() {
        let v1 = id_to_idx[&face.vertices[0]];
        let v2 = id_to_idx[&face.vertices[1]];
        let v3 = id_to_idx[&face.vertices[2]];
        xml.push_str(&format!(
            "        <triangle v1=\"{v1}\" v2=\"{v2}\" v3=\"{v3}\" />\n"
        ));
    }
    xml.push_str("      </triangles>\n");

    xml.push_str(MODEL_FOOTER);
    xml
}

const CONTENT_TYPES_XML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml" />
  <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml" />
</Types>
"#;

const RELS_XML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel" />
</Relationships>
"#;

const MODEL_HEADER: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<model unit="millimeter" xml:lang="en-US"
       xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">
  <resources>
    <object id="1" type="model">
      <mesh>
"#;

const MODEL_FOOTER: &str = r#"      </mesh>
    </object>
  </resources>
  <build>
    <item objectid="1" />
  </build>
</model>
"#;
