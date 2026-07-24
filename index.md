# MethScope

<div class="ms-wrap">
<p class="ms-tag">Analysis of DNA methylomes via Most Recurrent Methylation Patterns (MRMPs) — annotate, deconvolve, impute, and embed bulk, single-cell &amp; spatial data, uniquely robust on <b>sparse &amp; ultra-sparse</b> inputs.</p>
<div class="ms-badges"><span class="ms-badge"><b>C</b> CLI + <b>R</b> package</span><span class="ms-badge">CRAN <b>1.0.3</b></span><span class="ms-badge">R &ge; 4.0</span><span class="ms-badge">AGPL-3.0</span><span class="ms-badge">no-GPU</span></div>
<div class="ms-cta"><a class="external-link ms-btn primary" href="https://zhou-lab.github.io/methscope-cli/">⌨ Command line <span class="k">methscope-cli · C, no R</span></a><a class="ms-btn ghost" href="https://cran.r-project.org/package=MethScope">📦 R package</a><a class="ms-btn ghost" href="https://github.com/zhou-lab/methscope_data">🗂 Pretrained models</a></div>
<div class="ms-sech" id="applications">Applications <span class="sub">what MethScope does</span></div>
<div class="ms-grid">
<a class="ms-card" href="articles/MethScope-Tutorial.html"><span class="ic">🧬</span><span class="h">Sparse methylome annotation</span><span class="p">Predict cell type, sex, age and other traits from sparse methylomes, with confidence.</span><span class="go">predict →</span></a>
<a class="ms-card" href="articles/MethScope-Tutorial.html"><span class="ic">🩸</span><span class="h">Deconvolution</span><span class="p">Estimate cell-type proportions in mini-bulk / cfDNA mixtures (NNLS).</span><span class="go">deconv →</span></a>
<a class="ms-card" href="articles/MethScope-Tutorial.html"><span class="ic">🧩</span><span class="h">Imputation &amp; upscaling</span><span class="p">Reconstruct CpG-level methylation from ultra-sparse input.</span><span class="go">upscale →</span></a>
<a class="ms-card" href="articles/MethScope-Tutorial.html"><span class="ic">🗺️</span><span class="h">Clustering &amp; embedding</span><span class="p">MRMP representations for UMAP, clustering, and structure discovery.</span><span class="go">matrix →</span></a>
</div>
<div class="ms-sech">Two ways to run <span class="sub">same MRMP models, shell or R</span></div>
<div class="ms-two">
<div class="ms-panel"><div class="ph">⌨ Command line <span class="ms-tagp cli">methscope-cli · no R runtime</span></div><pre class="ms-pre dark"><span class="c"># install — conda package coming soon (no R runtime)</span>
conda install -c bioconda methscope-cli
<span class="c"># annotate · deconvolve · impute</span>
methscope <span class="p">predict</span> query.cg model.ubjx &gt; labels.tsv
methscope <span class="p">deconv</span>  mixture.cg panel.refx &gt; props.tsv
methscope <span class="p">upscale</span> -o out.cg model.updecx query.cg</pre><div class="ms-foot"><a class="external-link" href="https://zhou-lab.github.io/methscope-cli/">▸ runnable examples →</a></div></div>
<div class="ms-panel"><div class="ph">📦 R package <span class="ms-tagp r">interactive · plotting</span></div><pre class="ms-pre lite"><span class="c"># install (CRAN)</span>
install.packages("MethScope")
<span class="c"># annotate + visualize in R</span>
x &lt;- <a class="fn" href="reference/GenerateInput.html">GenerateInput</a>(cg, mrmp)
m &lt;- <a class="fn" href="reference/Zhou2025_HumanAtlas_P1000.html">Zhou2025_HumanAtlas_P1000</a>()
p &lt;- <a class="fn" href="reference/PredictCellType.html">PredictCellType</a>(m, x)
<a class="fn" href="reference/PlotUMAP.html">PlotUMAP</a>(x, p)</pre><div class="ms-foot"><a href="articles/MethScope-Tutorial.html">▸ runnable example →</a></div></div>
</div>
<div class="ms-sech">Documentation <span class="sub">go deeper</span></div>
<div class="ms-docs">
<a class="ms-doc" href="articles/MethScope-Input.html"><span class="t">Prepare input files →</span><span class="d">BED / ALLC / bigwig → .cg</span></a>
<a class="ms-doc" href="articles/MethScope-MRMP.html"><span class="t">Build &amp; interpret MRMPs →</span><span class="d">recurrent methylation patterns</span></a>
<a class="ms-doc" href="articles/pretrained-models.html"><span class="t">Work with pretrained models →</span><span class="d">bundle / unbundle &amp; formats</span></a>
<a class="ms-doc" href="reference/index.html"><span class="t">Function reference →</span><span class="d">all exported functions</span></a>
</div>
</div>
