import subprocess
import os

DOT_ARCH = """
digraph SHINE {
    /* ---- global style ---- */
    graph [fontname="Helvetica", rankdir=TB, nodesep=0.6, ranksep=0.8,
           bgcolor="transparent"];
    node  [fontname="Helvetica", shape=box, style="rounded,filled",
           fontsize=12, margin="0.25,0.12"];
    edge  [fontname="Helvetica", fontsize=10, color="#555555"];

    /* ---- top-level package ---- */
    shine [label="shine\\n(package root)", fillcolor="#4A90D9",
           fontcolor=white, fontsize=14, penwidth=2];

    /* ---- core modules ---- */
    core  [label="SHINE.py\\nCore extraction engine",
           fillcolor="#5DADE2", fontcolor=white];
    utils [label="shine_utils.py\\nUtility functions",
           fillcolor="#5DADE2", fontcolor=white];

    /* ---- sub-packages ---- */
    fem   [label="Find_Em_SHINE\\nLine-emitter pipeline",
           fillcolor="#27AE60", fontcolor=white, fontsize=13, penwidth=2];
    mim   [label="Make_Im_SHINE\\n2-D image generation",
           fillcolor="#E67E22", fontcolor=white, fontsize=13, penwidth=2];

    /* ---- Find_Em_SHINE modules ---- */
    ext   [label="extraction.py\\nStep 1 – SHINE extraction",
           fillcolor="#82E0AA"];
    cov   [label="covariance.py\\nStep 2 – Noise covariance",
           fillcolor="#82E0AA"];
    cat   [label="catalogue.py\\nStep 3 – Catalogue & cutouts",
           fillcolor="#82E0AA"];
    cli   [label="cli.py\\nCommand-line interface",
           fillcolor="#82E0AA"];

    /* ---- Make_Im_SHINE modules ---- */
    mmod  [label="Make_Im_SHINE.py\\nImage creation",
           fillcolor="#F0B27A"];

    /* ---- edges ---- */
    shine -> core  [label=" core"];
    shine -> utils [label=" utils"];
    shine -> fem   [label=" sub-pkg"];
    shine -> mim   [label=" sub-pkg"];

    fem -> ext [label=" step 1"];
    fem -> cov [label=" step 2"];
    fem -> cat [label=" step 3"];
    fem -> cli [label=" CLI"];

    mim -> mmod;

    /* ---- cross-module dependencies ---- */
    edge [style=dashed, color="#999999"];
    ext  -> core  [label="calls runextraction"];
    ext  -> utils [label="calls filter_cube,\\nclean_clube"];
    cat  -> cov   [label="uses covariance\\ncorrection"];
}
"""

DOT_PIPELINE = """
digraph pipeline {
    graph [fontname="Helvetica", rankdir=LR, nodesep=0.5, ranksep=1.0,
           bgcolor="transparent"];
    node  [fontname="Helvetica", shape=box, style="rounded,filled",
           fontsize=11, margin="0.2,0.1"];
    edge  [fontname="Helvetica", fontsize=9];

    /* ---- optional pre-processing ---- */
    pre  [label="Pre-processing\\n(optional)\\nclean_clube\\nfilter_cube",
          fillcolor="#D5F5E3", style="rounded,filled,dashed"];

    /* ---- pipeline steps ---- */
    s1   [label="Step 1\\nextraction.py\\n──────────\\nSegmentation map\\nFiltered cubes\\nRaw catalogue",
          fillcolor="#82E0AA"];
    s2   [label="Step 2\\ncovariance.py\\n──────────\\nCovariance .npz\\nPolynomial fit\\nDiagnostic plot",
          fillcolor="#82E0AA"];
    s3   [label="Step 3\\ncatalogue.py\\n──────────\\nS/N-corrected cat.\\nConfidence classes\\nImage cutouts\\n1-D spectra",
          fillcolor="#82E0AA"];

    pre -> s1 [label="cubes"];
    s1  -> s2 [label="filtered cubes\\n+ labels"];
    s1  -> s3 [label="raw catalogue\\n+ seg. map"];
    s2  -> s3 [label="covariance\\npolynomial"];
}
"""

def generate_svgs():
    docs_dir = "/data/dtornotti/SHINE/docs"
    static_dir = os.path.join(docs_dir, "_static")
    os.makedirs(static_dir, exist_ok=True)
    
    # Save files
    with open("/tmp/shine_arch.dot", "w") as f:
        f.write(DOT_ARCH)
    with open("/tmp/pipeline_flow.dot", "w") as f:
        f.write(DOT_PIPELINE)
        
    # Compile
    subprocess.run(["dot", "-Tsvg", "/tmp/shine_arch.dot", "-o", os.path.join(static_dir, "architecture_tree.svg")], check=True)
    subprocess.run(["dot", "-Tsvg", "/tmp/pipeline_flow.dot", "-o", os.path.join(static_dir, "pipeline_flow.svg")], check=True)
    print("SVGs generated successfully.")

if __name__ == "__main__":
    generate_svgs()
