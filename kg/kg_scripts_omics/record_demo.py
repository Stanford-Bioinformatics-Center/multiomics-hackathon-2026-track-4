"""Record the captioned explorer demo (Playwright + Chromium): python3 scripts/record_demo.py, then
ffmpeg -i demo_video/*.webm -filter:v "setpts=PTS/1.55,fps=30" -c:v libx264 -crf 22 -pix_fmt yuv420p demo.mp4
"""
"""Scripted, captioned walkthrough of the Exercise KG Explorer, recorded with Playwright."""
import asyncio, subprocess, time, os, glob, shutil
from playwright.async_api import async_playwright

EXP = os.environ.get("EXPLORER_DIR", "explorer")
OUT = os.environ.get("VIDEO_DIR", "demo_video")
W, H = 1440, 900
shutil.rmtree(OUT, ignore_errors=True); os.makedirs(OUT)
srv = subprocess.Popen(["python3", "-m", "http.server", "8777"], cwd=EXP, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
time.sleep(1)

OVERLAY = r"""
<style>
#demo-cap{position:fixed;left:50%;bottom:26px;transform:translateX(-50%);max-width:760px;z-index:99999;
  background:rgba(12,20,18,.88);color:#fff;font:500 17px/1.4 -apple-system,"Segoe UI",Roboto,sans-serif;padding:12px 18px;border-radius:10px;
  box-shadow:0 8px 30px rgba(0,0,0,.25);transition:opacity .35s;opacity:0;pointer-events:none}
#demo-cap b{color:#7fe0cf;font-weight:600}
#demo-cap small{display:block;font-size:13px;color:#bcd;margin-top:3px;font-weight:400}
#demo-cur{position:fixed;left:0;top:0;width:22px;height:22px;z-index:100000;pointer-events:none;transition:transform .55s cubic-bezier(.3,.7,.3,1)}
#demo-cur svg{filter:drop-shadow(0 1px 2px rgba(0,0,0,.4))}
#demo-ring{position:fixed;width:34px;height:34px;margin:-17px 0 0 -17px;border:3px solid #0E6E62;border-radius:50%;z-index:99998;pointer-events:none;opacity:0}
#demo-ring.go{animation:ring .5s ease-out}
@keyframes ring{0%{opacity:.9;transform:scale(.3)}100%{opacity:0;transform:scale(1.4)}}
#demo-title{position:fixed;inset:0;z-index:100001;background:#0F1513;color:#E3ECE8;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:14px;
  font-family:-apple-system,"Segoe UI",Roboto,sans-serif;transition:opacity .6s;opacity:0;pointer-events:none;text-align:center;padding:40px}
#demo-title h1{font-size:44px;margin:0;font-weight:650;letter-spacing:.01em}
#demo-title p{font-size:20px;margin:0;color:#94A49F;max-width:900px;line-height:1.5}
#demo-title .k{color:#4FB8A8}
</style>
<div id="demo-cap"></div><div id="demo-ring"></div>
<div id="demo-cur"><svg width="22" height="22" viewBox="0 0 22 22"><path d="M2 2 L2 18 L7 13.5 L10.5 20.5 L13.2 19.2 L9.8 12.4 L16.5 12.4 Z" fill="#fff" stroke="#111" stroke-width="1.4" stroke-linejoin="round"/></svg></div>
<div id="demo-title"></div>
"""


async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        ctx = await b.new_context(viewport={"width": W, "height": H}, record_video_dir=OUT, record_video_size={"width": W, "height": H})
        pg = await ctx.new_page()
        await pg.route("**/cytoscape.min.js", lambda r: r.fulfill(path=f"{EXP}/vendor/cytoscape.min.js", content_type="application/javascript"))
        await pg.route("**/fonts.googleapis.com/**", lambda r: r.abort())
        html = open(f"{EXP}/index.html", encoding="utf-8").read()
        k = html.rindex("})();")
        html = '<!doctype html><meta charset="utf-8">' + html[:k] + "window.__T={focusOn, get D(){return D}, get cy(){return cy}};" + html[k:] + OVERLAY
        await pg.route("**/demo.html", lambda r: r.fulfill(body=html, content_type="text/html; charset=utf-8"))

        cur = [W // 2, H // 2]

        async def title(h1, sub, secs):
            await pg.evaluate("([h,s])=>{const t=document.getElementById('demo-title'); t.innerHTML=`<h1>${h}</h1><p>${s}</p>`; t.style.opacity=1}", [h1, sub])
            await pg.wait_for_timeout(int(secs * 1000))
            await pg.evaluate("document.getElementById('demo-title').style.opacity=0")
            await pg.wait_for_timeout(600)

        async def cap(text, secs=0.0):
            await pg.evaluate("t=>{const c=document.getElementById('demo-cap'); if(!t){c.style.opacity=0;return;} c.innerHTML=t; c.style.opacity=1}", text)
            if secs: await pg.wait_for_timeout(int(secs * 1000))

        async def move(x, y, pause=0.65):
            cur[0], cur[1] = x, y
            await pg.evaluate("([x,y])=>{document.getElementById('demo-cur').style.transform=`translate(${x}px,${y}px)`}", [x, y])
            await pg.mouse.move(x, y)
            await pg.wait_for_timeout(int(pause * 1000))

        async def ring(x, y):
            await pg.evaluate("([x,y])=>{const r=document.getElementById('demo-ring'); r.style.left=x+'px'; r.style.top=y+'px'; r.classList.remove('go'); void r.offsetWidth; r.classList.add('go')}", [x, y])

        async def click(sel, after=0.9, nth=0):
            loc = pg.locator(sel).nth(nth)
            try: await loc.scroll_into_view_if_needed(timeout=4000)
            except Exception as e:
                print("MISSING click", sel, nth); await pg.screenshot(path=f"{OUT}/miss_{len(sel)}.png"); return
            bb = await loc.bounding_box()
            x, y = bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2
            await move(x, y)
            await ring(x, y)
            await pg.mouse.click(x, y)
            await pg.wait_for_timeout(int(after * 1000))

        async def chip(container, text, after=0.9):
            await click(f"#{container} .chip:text-is('{text}')", after)

        async def slider(sel, value, after=1.2):
            loc = pg.locator(sel)
            await loc.scroll_into_view_if_needed(timeout=4000)
            bb = await loc.bounding_box()
            mn, mx = [float(await loc.get_attribute(a)) for a in ("min", "max")]
            x = bb["x"] + 8 + (bb["width"] - 16) * (value - mn) / (mx - mn)
            y = bb["y"] + bb["height"] / 2
            await move(x, y, .4)
            await ring(x, y)
            await pg.evaluate("([s,v])=>{const e=document.querySelector(s); e.value=v; e.dispatchEvent(new Event('input',{bubbles:true}))}", [sel, value])
            await pg.wait_for_timeout(int(after * 1000))

        async def node(js_pred, after=1.4):
            """click a graph node chosen by a JS predicate over node data d"""
            pos = await pg.evaluate("""pred=>{const f=new Function('d','D','return ('+pred+')'); const D=__T.D, cy=__T.cy;
                 const n=cy.nodes().filter(n=>f(n.data(),D))[0]; if(!n) return null; const r=document.getElementById('cy').getBoundingClientRect();
                 const p=n.renderedPosition(); return [r.left+p.x, r.top+p.y];}""", js_pred)
            if not pos:
                print("node not found:", js_pred); return False
            await move(*pos)
            await ring(*pos)
            await pg.mouse.click(*pos)
            await pg.wait_for_timeout(int(after * 1000))
            return True

        async def edge(js_pred, after=1.4):
            pos = await pg.evaluate("""pred=>{const f=new Function('d','D','return ('+pred+')'); const D=__T.D, cy=__T.cy;
                 const e=cy.edges().filter(e=>f(e.data(),D))[0]; if(!e) return null; const r=document.getElementById('cy').getBoundingClientRect();
                 const p=e.renderedMidpoint(); return [r.left+p.x, r.top+p.y];}""", js_pred)
            if not pos:
                print("edge not found:", js_pred); return False
            await move(*pos); await ring(*pos); await pg.mouse.click(*pos)
            await pg.wait_for_timeout(int(after * 1000))
            return True

        async def search(text, pick, after=2.2):
            await click("#q", .3)
            await pg.fill("#q", "")
            await pg.type("#q", text, delay=90)
            await pg.wait_for_timeout(600)
            await click(f"#sugg button:has(span:text-is('{pick}'))", after)

        async def hover(sel, secs=1.2, nth=0):
            loc = pg.locator(sel).nth(nth)
            try: await loc.scroll_into_view_if_needed(timeout=4000)
            except Exception:
                print("MISSING hover", sel, nth); await pg.screenshot(path=f"{OUT}/missh_{nth}.png"); return
            bb = await loc.bounding_box()
            await move(bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2, secs)

        async def scroll_pane(sel, px, secs=1.2):
            await pg.evaluate("([s,p])=>{const e=document.querySelector(s); e.scrollBy({top:p,behavior:'smooth'})}", [sel, px])
            await pg.wait_for_timeout(int(secs * 1000))

        # ------------------------------------------------------------------ load
        await pg.goto("http://localhost:8777/demo.html")
        await pg.evaluate("document.getElementById('demo-title').style.opacity=1")
        await pg.wait_for_function("window.__T && __T.D && __T.D.genes && document.querySelectorAll('#hyp-list .card').length>0", timeout=60000)
        await pg.wait_for_timeout(1500)

        # ------------------------------------------------------------------ 0. title
        await title("Exercise KG Explorer",
                    "MoTrPAC multi-omics as a knowledge graph<br><span class='k'>Rat</span> 8-week endurance training · 20 organs &nbsp;|&nbsp; "
                    "<span class='k'>Human</span> one endurance or resistance bout · muscle, adipose, blood<br>"
                    "Significant responses only · pathways · protein interactions · phenotypes · a citation on every link", 6)

        # ------------------------------------------------------------------ 1. overview
        await cap("<b>Pathway view.</b> A pathway, its responding genes and molecules, and the organs where they change.<small>Red = up, blue = down; edge width = strength (percentile of −log10 adjusted p within assay).</small>", 4.5)
        await hover("#legend", 2.0)
        await cap("<b>Each node keeps its 5 strongest edges</b> by default, so the view stays readable.", 3)

        # ------------------------------------------------------------------ 2. filters
        await cap("<b>Filters.</b> Species…", 0.5)
        await chip("f-species", "Rat", 1.6); await chip("f-species", "Rat", 1.0)
        await cap("<b>Filters.</b> …timepoint, sex, direction…", 0.3)
        await chip("f-tp", "1 wk", .8); await chip("f-tp", "2 wk", 1.0); await click("#tp-all", 1.0)
        await chip("f-dir", "Down", 1.2); await chip("f-dir", "Down", .8)
        await cap("<b>Filters.</b> …omics layer and interaction type…", 0.3)
        await chip("f-layer", "Transcript", 1.2); await chip("f-layer", "Transcript", .8)
        await chip("f-ix", "Protein–protein", 1.2); await chip("f-ix", "Protein–protein", .8)
        await cap("<b>Minimum strength</b> keeps only the strongest responses; <b>connect to</b> organs or exercise groups.", 0.3)
        await slider("#f-z", 60, 1.6); await slider("#f-z", 0, .8)
        await chip("f-target", "Exercise groups", 1.8); await chip("f-target", "Organs", 1.0)
        await cap("<b>Max edges per node</b> (5 by default) and how many molecules to draw.", 0.3)
        await slider("#f-deg", 12, 1.6); await slider("#f-deg", 5, .8)
        await cap("<b>Layouts:</b> force-directed or rings (pathway → genes → molecules → organs).", 0.3)
        await click("#btn-layout", 2.0); await click("#btn-layout", 1.2)

        # ------------------------------------------------------------------ 3. muscle: SLC7A5
        await cap("<b>Muscle · SLC7A5</b> (LAT1, the leucine transporter). Gene view, muscle preset.", 0.3)
        await click("#tab-gene", 1.0)
        await click("[data-preset=muscle]", .6)
        await search("SLC7A5", "SLC7A5", 2.4)
        await cap("<b>SLC7A5 in human muscle:</b> the transcript rises 3.5–4 h after both endurance and resistance exercise, then dips below control 24 h after resistance.", 3.5)
        await node("d.kind==='feature' && D.features.id[d.ref].startsWith('human:transcript')", 1.0)
        await cap("<b>Details:</b> every significant response (organ, group, timepoint, adjusted p, z) plus database links for the molecule and its gene.", 3.2)
        await hover("#pane-details .xrefs .xref", 1.4)
        await hover("#pane-details .xrefs .xref", 1.0, nth=2)
        await scroll_pane("#pane-details", 260, 1.6)
        await cap("<b>Sources</b> for everything in the panel, with DOIs.", 0.3)
        await click("#pane-details details.reflist summary", 1.8)
        ok = await edge("d.kind==='ppi'", 1.0)
        if ok: await cap("<b>Edges are clickable too:</b> a physical interaction from STRING v12 / Hetionet with its confidence score.", 3)
        await click("#ptab-hyp", .6)
        await cap("<b>Hypotheses tab.</b> Green = what the graph shows; amber = new hypotheses. Every sentence carries its own citation.", 4)
        await scroll_pane("#pane-hyp", 380, 2.2)
        await scroll_pane("#pane-hyp", -380, .8)

        # ------------------------------------------------------------------ 4. rat organs for the same gene
        await cap("<b>Same gene, all organs:</b> rat Slc7a5 rises with training in brown fat, blood, colon and spleen.", 0.3)
        await click("[data-preset=all]", 3.2)
        await node("d.kind==='feature' && D.features.id[d.ref].startsWith('rat:')", 2.8)

        # ------------------------------------------------------------------ 5. adipose
        await cap("<b>Adipose.</b> NR4A1 is the strongest human adipose transcript response to one bout.", 0.3)
        await click("[data-preset=adipose]", .6)
        await search("NR4A1", "NR4A1", 2.6)
        await node("d.kind==='feature' && D.features.id[d.ref].startsWith('human:')", 2.6)

        # ------------------------------------------------------------------ 6. blood
        await cap("<b>Blood · metabolites.</b> Metabolites link only to organs, other metabolites and proteins.", 0.3)
        await click("#tab-metabolite", .6)
        await click("[data-preset=blood]", .6)
        await search("lactic", "Lactic acid", 2.6)
        await node("d.kind==='feature' && D.features.type[d.ref]===6 && D.features.name[d.ref]==='Lactic acid'", 1.2)
        await cap("<b>Lactic acid</b> climbs in human blood during and just after exercise; linked to KEGG and PubChem.", 3.0)
        await click("#ptab-hyp", .5)
        await cap("<b>Co-regulated proteins</b> (dashed) are hypotheses about producers or consumers, cited as such.", 2.4)

        # ------------------------------------------------------------------ 7. every rat organ
        await cap("<b>All 23 organs.</b> Lactic acid across every rat organ and human tissue with data.", 0.3)
        await click("[data-preset=all]", 3.4)
        await node("d.kind==='target' && d.tkind==='tissue' && d.label==='Rat Heart'", 2.2) or await node("d.kind==='target'", 2.2)

        # ------------------------------------------------------------------ 8. pathway + kinase hypotheses
        await cap("<b>Pathways with kinases:</b> AMPK signalling, phosphosites and their upstream kinases (dotted).", 0.3)
        await click("#tab-pathway", .6)
        await click("[data-preset=muscle]", .5)
        await search("ampk", "ACTIVATION OF AMPK DOWNSTREAM OF NMDARS", 3.0)
        await node("d.kind==='gene' && d.kinase===1", 1.4) or await node("d.kind==='pathway'", 1.4)
        await click("#ptab-hyp", .5)
        await cap("<b>Hypothesis cards</b> infer kinase activity from substrate sites and name a test, with literature citations.", 3.4)
        await scroll_pane("#pane-hyp", 300, 1.6)
        await click("#hyp-list .ent", 2.6)

        # ------------------------------------------------------------------ 9. phenotypes
        await cap("<b>Phenotype view:</b> rat VO2max, body composition and muscle mass; human clinical chemistry.", 0.3)
        await click("#tab-phenotype", 3.0)
        await node("d.kind==='pheno' && d.label==='Insulin'", 1.6) or await node("d.kind==='pheno'", 1.6)
        await cap("<b>Clinical markers</b> link to NCBI Gene / UniProt / KEGG; rat phenotypes use Welch t-tests vs sedentary controls.", 3.2)
        await node("d.kind==='pheno' && d.label==='VO2max change'", 2.4)

        # ------------------------------------------------------------------ 9b. themes and the disease layer
        async def theme(val, after):
            loc = pg.locator("#theme"); bb = await loc.bounding_box()
            x, y = bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2
            await move(x, y); await ring(x, y)
            await pg.select_option("#theme", val); await pg.wait_for_timeout(int(after * 1000))
        # --- cancer (theme) with disease details
        await cap("<b>Cancer.</b> Themes load a disease with the organs that matter. Colon cancer: its genes are over-represented among exercise responders in rat colon.", 0.3)
        await theme("cancer", 3.2)
        await node("d.kind==='disease'", 1.4)
        await cap("<b>Disease details:</b> Disease Ontology link, over-representation per organ (hypergeometric, BH) and the drugs that treat it.<small>Disease layer from Hetionet (Disease Ontology, DrugBank).</small>", 3.2)
        await edge("d.kind==='dcon' && d.opp===1", 1.4)
        await cap("<b>Direction test:</b> exercise moves colon-cancer genes <b>against</b> the disease in human blood (green) but <b>with</b> it in rat colon (amber).", 3.4)
        await click("#ptab-hyp", .5)
        await cap("<b>Hypotheses:</b> opposing a disease signature is one route to exercise as therapy; drugs whose targets exercise also moves are flagged. All cited.", 3.4)

        # --- brain (theme)
        await cap("<b>Brain.</b> Schizophrenia genes are over-represented among exercise responders in rat cortex and hippocampus.", 0.3)
        await theme("brain", 3.2)
        pos = await pg.evaluate("""()=>{const cy=__T.cy; const e=cy.edges().filter(e=>e.data('kind')==='greg' && /Cortex|Hippocampus/.test(e.target().data('label')))[0];
             if(!e) return null; const n=e.source(); const r=document.getElementById('cy').getBoundingClientRect(); const p=n.renderedPosition(); return [r.left+p.x, r.top+p.y, n.data('label'), e.target().data('label')];}""")
        if pos:
            await move(pos[0], pos[1]); await ring(pos[0], pos[1]); await pg.mouse.click(pos[0], pos[1]); await pg.wait_for_timeout(900)
            await cap(f"<b>{pos[2]}</b> is linked to schizophrenia and responds to training in rat {pos[3].replace('Rat ','').lower()}. Drugs for schizophrenia that bind responding genes are shown in blue.", 3.2)

        # --- gut (gene view)
        await cap("<b>Gut.</b> PYY, a gut hormone that signals fullness, rises in the trained rat small intestine in weeks 1–2 (and in brown fat and adrenal).", 0.3)
        await click("#tab-gene", .5)
        await click("[data-preset=all]", .5)
        await search("Pyy", "Pyy", 2.6)
        await node("d.kind==='feature' && D.features.id[d.ref]==='rat:IMMUNO:PYY'", 2.6)
        await cap("<b>Gut and appetite:</b> in the same rats plasma leptin falls with training, while leptin in the small intestine rises.", 0.3)
        await search("Lep", "Lep", 2.4)
        await node("d.kind==='feature' && D.features.id[d.ref]==='rat:IMMUNO:LEPTIN'", 2.6)

        # ------------------------------------------------------------------ 10. Claude
        await cap("<b>Ask Claude</b> and <b>Generate insights</b> query the graph with tools and cite every sentence.<small>These run when the explorer is opened inside Claude; this recording is the standalone Docker build.</small>", 0.3)
        await click("#ptab-ask", 4.5)
        await click("#ptab-hyp", .6)
        await hover("#btn-insights", 2.0)
        await cap("", 0)

        # ------------------------------------------------------------------ end
        await title("Exercise KG",
                    "77,123 nodes · 1.50 M relationships · 53 cited sources · diseases &amp; drugs<br>"
                    "<span class='k'>bash deploy.sh</span> → Neo4j on :7474 and this explorer on :8080<br>"
                    "Methods & caveats: docs/METHODS.md · adjacency matrices: kg/adjacency/", 5)
        await ctx.close()
        await b.close()


asyncio.run(main())
srv.terminate()
print(glob.glob(f"{OUT}/*.webm"))
