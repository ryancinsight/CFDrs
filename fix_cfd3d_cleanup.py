"""
Targeted fix for cfd-3d corruption from previous regex scripts.
Strategy:
1. Remove `: cfd_mesh::domain::core::Scalar` from INSIDE type-parameter applications
   (e.g. Option<T: Scalar> -> Option<T>, Vec<T: Scalar> -> Vec<T>, etc.)
2. At impl/struct/enum DECLARATION lines: ensure T has the Scalar bound when T is the
   parameter being declared (not applied).
"""
import os
import re

SCALAR_BOUND = "cfd_mesh::domain::core::Scalar"
SCALAR_TAG   = f": {SCALAR_BOUND}"

def strip_bad_tag(content: str) -> str:
    """
    Remove `: cfd_mesh::domain::core::Scalar` that appears INSIDE 
    a type-parameter application like Option<T: Scalar>, Vec<T: Scalar>, etc.
    These are E0229 ("associated item constraint not allowed here") errors.
    
    We remove any occurrence of `: cfd_mesh::domain::core::Scalar` that is
    followed by `>` (closing a type application) or `,` (in a generic list).
    """
    # Remove tag inside type applications: <T: cfd_mesh::domain::core::Scalar>
    # but NOT at top-level struct/impl/enum declarations.
    # Pattern: any word<T: SCALAR_BOUND[, ...]>
    # We'll strip the tag in field type positions by replacing:
    #   SomeName<T: cfd_mesh::domain::core::Scalar>  ->  SomeName<T>
    #   SomeName<T: cfd_mesh::domain::core::Scalar, U>  ->  SomeName<T, U>
    pattern = re.compile(
        r'(?<![a-z _(>])'        # not preceded by lowercase / normal decl context
        r'((?:(?:[\w:]+)\s*<[^<>]*?)'  # something like SomeName<...
        r'T): cfd_mesh::domain::core::Scalar'  # T: SCALAR
        r'(?=[^a-zA-Z])',              # followed by non-alpha (>, ,, etc)
    )
    # Simpler: just remove ": cfd_mesh::domain::core::Scalar" whenever it appears  
    # inside angle brackets that are part of a type application (heuristic: preceded
    # by a generic type name token, not "struct", "enum", "impl", "pub", "fn").
    
    # Replace `<T: cfd_mesh::domain::core::Scalar>` → `<T>`  ONLY when preceded by
    # an identifier that's a type name (not a keyword like `struct`/`impl`/`enum`).
    # We detect "type application" vs "declaration" by whether the angle bracket is
    # preceded by an IDENTIFIER (type usage) vs a keyword (declaration).
    
    # Step 1: handle the common patterns from the corruption
    # Pattern: Identifier<T: SCALAR_BOUND> → Identifier<T>
    content = re.sub(
        r'(\b(?:Option|Vec|DMatrix|DVector|Matrix|Point3|Vector3|HashMap|'
        r'BoundaryCondition|ConstantPropertyFluid|SolverConfig|IndexedMesh|'
        r'BifurcationGeometry3D|BifurcationConfig3D|BifurcationMesh|'
        r'MeshRefinementConfig|ConicalTransition|ForcingMethod|'
        r'SerpentineConfig3D|SerpentineSolution3D|TrifurcationConfig3D|'
        r'TrifurcationSolution3D|VenturiConfig3D|VenturiSolution3D|'
        r'StokesFlowProblem|ConstantFluidConfig|NavierStokesFluid|'
        r'SolveResult|NavierStokesVOF|BranchingMeshBuilder|'
        r'SerpentineMeshBuilder|VenturiMeshBuilder|TrifurcationMeshBuilder|'
        r'VolumeMesh|HierarchicalMesh|TriFace|TetCell|'
        r'Cell|Face|Vertex|BoxLayout|[A-Z][a-zA-Z0-9_]*)'
        r')<T: cfd_mesh::domain::core::Scalar>',
        r'\1<T>',
        content
    )
    
    # Step 2: Box<dyn Trait<T: SCALAR>> → Box<dyn Trait<T>>
    content = re.sub(
        r'(Box<dyn \w+)<T: cfd_mesh::domain::core::Scalar>',
        r'\1<T>',
        content
    )
    
    # Step 3: impl<T: SCALAR + REST> TypeName<T: SCALAR> → impl<T: SCALAR + REST> TypeName<T>
    # Remove the bogus `: SCALAR` from the TYPE APPLICATION (after struct/impl name)
    content = re.sub(
        r'(impl<[^>]+>)\s+([A-Za-z0-9_:]+)<T: cfd_mesh::domain::core::Scalar>',
        r'\1 \2<T>',
        content
    )
    
    # Step 4: struct/enum/pub struct XYZ<T: cfd_mesh::domain::core::Scalar> stays OK
    # (these are declaration sites and the bound is correct there)
    # But if the struct was NOT originally generic, the script added spurious bounds.
    # We can detect "struct Foo<T: SCALAR + RealField + Copy>" which is desired.
    
    # Step 5: fn foo() -> SomeType<T: SCALAR> → fn foo() -> SomeType<T>
    content = re.sub(
        r'(->\s*[A-Za-z0-9_:<>]*?)<T: cfd_mesh::domain::core::Scalar>',
        r'\1<T>',
        content
    )
    
    # Step 6: param: SomeType<T: SCALAR> in function args
    content = re.sub(
        r'(\w: [A-Za-z0-9_:<>]*?)<T: cfd_mesh::domain::core::Scalar>',
        r'\1<T>',
        content
    )
    
    return content


def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        original = f.read()
    
    fixed = strip_bad_tag(original)
    
    if fixed != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(fixed)
        print(f"  Fixed: {filepath}")


def main():
    root = 'crates/cfd-3d/src'
    count = 0
    for dirpath, _, files in os.walk(root):
        for fname in files:
            if fname.endswith('.rs'):
                process_file(os.path.join(dirpath, fname))
                count += 1
    print(f"\nProcessed {count} files.")


if __name__ == '__main__':
    main()
