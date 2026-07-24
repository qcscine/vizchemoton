"""
Enric Petrus, December 2024. Added SCINE helper function to link with the
amk-tools generation of HTML files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link
amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
"""

# Standard library imports
from collections import Counter
from operator import itemgetter
import copy

# Third-party library imports
import numpy as np
import bokeh.plotting
import bokeh.models as bkm
import RXVisualizer as arxviz
import RXReader as arx
import networkx as nx
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# HTML definitions

# SMILES editor full, minified iframe
iframe_editor = """<iframe width="100%" height="1000px" frameBorder="0" srcdoc="<!DOCTYPE html><html lang=&quot;en&quot;><head><meta charset=&quot;UTF-8&quot;><meta name=&quot;viewport&quot; content=&quot;width=device-width, initial-scale=1.0&quot;><title>2D structure ⇄ SMILES/InChIKey</title><link rel=&quot;stylesheet&quot; type=&quot;text/css&quot; href=&quot;https://cdn.jsdelivr.net/npm/kekule/dist/themes/default/kekule.css&quot; /><script src=&quot;https://cdn.jsdelivr.net/npm/kekule/dist/kekule.min.js?modules=chemWidget,algorithm,io&quot;></script><script src=&quot;https://unpkg.com/@rdkit/rdkit/Dist/RDKit_minimal.js&quot;></script><style> body { font-family: -apple-system, BlinkMacSystemFont, &quot;Segoe UI&quot;, Roboto, sans-serif; margin: 20px; background-color: #ffffff; } .container { max-width: 1000px; margin: 0 auto; } .page-title { text-align: center; color: #1f2937; margin-bottom: 25px; } .info-box { padding: 12px 15px; margin-bottom: 12px; border-radius: 4px; font-size: 0.95em; } .instructions-box { background-color: #f3f4f6; border-left: 4px solid #9ca3af; color: #374151; cursor: pointer; } .instructions-box summary { font-weight: bold; outline: none; user-select: none; } .instructions-box ul { margin: 10px 0 0 0; padding-left: 20px; } .privacy-note { background-color: #eff6ff; border-left: 4px solid #3b82f6; color: #1e3a8a; margin-bottom: 20px; } #composer { width: 100%; height: 380px; border: 1px solid #ccc; border-radius: 4px; background-color: #fff; } .btn-container { margin-top: 25px; margin-bottom: 25px; display: flex; gap: 20px; justify-content: center; align-items: center; } .calc-btn { color: white; border: none; padding: 14px 20px; font-size: 1.1em; font-weight: bold; border-radius: 4px; cursor: pointer; transition: background 0.2s, transform 0.1s; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); flex: 2; max-width: 320px; text-align: center; box-sizing: border-box; } .btn-clear { color: #4b5563; background-color: #ffffff; border: 2px solid #4b5563; padding: 12px 20px; font-size: 1.1em; font-weight: bold; border-radius: 4px; cursor: pointer; transition: background 0.2s, transform 0.1s; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); flex: 1; max-width: 160px; text-align: center; box-sizing: border-box; white-space: nowrap; } .calc-btn:active, .btn-clear:active { transform: scale(0.98); } .btn-2d-smiles { background-color: #7c3aed; } .btn-2d-smiles:hover { background-color: #6d28d9; } .btn-smiles-2d { background-color: #ea580c; } .btn-smiles-2d:hover { background-color: #c2410c; } .btn-clear:hover { background-color: #f3f4f6; } .calc-btn:disabled, .btn-clear:disabled { background-color: #cccccc !important; color: #666666 !important; border-color: #cccccc !important; cursor: not-allowed; transform: none !important; box-shadow: none; } .status-message { text-align: center; font-weight: 600; font-size: 1.05em; min-height: 24px; margin-bottom: 0px; } .status-success { color: #16a34a; } .status-error { color: #dc2626; } .output-panel { padding: 5px 15px 15px 15px; background-color: #fff; border: 1px solid #ddd; border-radius: 4px; } .output-group { margin-bottom: 15px; } .output-group:last-child { margin-bottom: 5px; } .output-row { display: flex; gap: 10px; align-items: stretch; margin-top: 5px; } .key-display { font-family: monospace; background-color: #f1f1f1; padding: 10px; border: 1px solid #e5e7eb; border-radius: 4px; font-size: 1.1em; word-break: break-all; min-height: 22px; flex: 1; display: flex; align-items: center; } textarea.key-display { resize: vertical; display: block; margin: 0; width: auto; } .copy-btn { background-color: #f3f4f6; color: #4b5563; border: 1px solid #d1d5db; border-radius: 4px; padding: 0 15px; font-weight: 600; cursor: pointer; transition: all 0.2s; display: flex; align-items: center; justify-content: center; min-width: 85px; user-select: none; } .copy-btn:hover { background-color: #e5e7eb; color: #1f2937; } .copy-btn.copied { background-color: #dcfce7; color: #16a34a; border-color: #bbf7d0; } .info-link { color: #007bff; text-decoration: none; } .info-link:hover { text-decoration: underline; } </style></head><body><div class=&quot;container&quot;><div id=&quot;composer&quot; data-widget=&quot;Kekule.Editor.Composer&quot; data-enable-style-toolbar=&quot;false&quot; data-common-tool-buttons=&quot;['undo', 'redo']&quot; data-chem-tool-buttons=&quot;['manipulate', 'erase', 'bond', 'atom', 'ring', 'charge']&quot;></div><div class=&quot;btn-container&quot;><button id=&quot;generate-btn&quot; class=&quot;calc-btn btn-2d-smiles&quot; disabled>Loading Core Engines...</button><button id=&quot;clear-btn&quot; class=&quot;btn-clear&quot; disabled>✕ Reset</button><button id=&quot;load-smiles-btn&quot; class=&quot;calc-btn btn-smiles-2d&quot; disabled>Loading Core Engines...</button></div><div id=&quot;status-msg&quot; class=&quot;status-message&quot;></div><div class=&quot;output-panel&quot;><div class=&quot;output-group&quot;><p style=&quot;margin: 10px 0 0 0;&quot;><strong>SMILES (Reference / Input):</strong></p><div class=&quot;output-row&quot;><textarea id=&quot;smiles-output&quot; class=&quot;key-display&quot; rows=&quot;2&quot; placeholder=&quot;Paste or generate SMILES here...&quot;></textarea><button class=&quot;copy-btn&quot; onclick=&quot;copyData('smiles-output', this)&quot;>📋 Copy</button></div></div><div class=&quot;output-group&quot;><p style=&quot;margin: 5px 0 0 0;&quot;><strong>InChIKey:</strong></p><div class=&quot;output-row&quot;><div id=&quot;inchikey-output&quot; class=&quot;key-display&quot;>Ready.</div><button class=&quot;copy-btn&quot; onclick=&quot;copyData('inchikey-output', this)&quot;>📋 Copy</button></div></div></div><hr style=&quot;border: 0; border-top: 1px solid #ddd; margin-top: 40px; margin-bottom: 20px;&quot;></div><script> let rdkitModule = null; let composerApp = null; const generateBtn = document.getElementById('generate-btn'); const loadSmilesBtn = document.getElementById('load-smiles-btn'); const clearBtn = document.getElementById('clear-btn'); const inchikeyOutput = document.getElementById('inchikey-output'); const smilesOutput = document.getElementById('smiles-output'); const statusMsg = document.getElementById('status-msg'); function showStatus(text, isSuccess) { statusMsg.innerText = text; statusMsg.className = &quot;status-message &quot; + (isSuccess ? &quot;status-success&quot; : &quot;status-error&quot;); } function copyData(elementId, buttonElement) { const target = document.getElementById(elementId); let textToCopy = target.tagName === &quot;TEXTAREA&quot; ? target.value : target.innerText; if (!textToCopy || textToCopy === &quot;Ready.&quot; || textToCopy.startsWith(&quot;Error:&quot;) || textToCopy.startsWith(&quot;Conversion failed:&quot;)) { return; } navigator.clipboard.writeText(textToCopy).then(() => { buttonElement.innerText = &quot;Copied! ✓&quot;; buttonElement.classList.add('copied'); setTimeout(() => { buttonElement.innerText = &quot;📋 Copy&quot;; buttonElement.classList.remove('copied'); }, 1500); }).catch(err => { console.error(&quot;Clipboard copy operation failed: &quot;, err); }); } function clearWholeApp() { if (composerApp) { composerApp.newDoc(); } inchikeyOutput.innerText = &quot;Ready.&quot;; smilesOutput.value = &quot;&quot;; statusMsg.innerText = &quot;&quot;; statusMsg.className = &quot;status-message&quot;; } window.initRDKitModule().then(function(instance) { rdkitModule = instance; checkInitializationComplete(); }).catch(err => { inchikeyOutput.innerText = &quot;Fatal Error: Unable to load RDKit WebAssembly core.&quot;; showStatus(&quot;Fatal Error: Unable to load RDKit core.&quot;, false); console.error(err); }); Kekule.X.domReady(() => { composerApp = Kekule.Widget.getWidgetById('composer'); checkInitializationComplete(); }); function checkInitializationComplete() { if (rdkitModule &amp;&amp; composerApp) { generateBtn.innerText = &quot;↓ SMILES&quot;; generateBtn.disabled = false; loadSmilesBtn.innerText = &quot;↑ 2D&quot;; loadSmilesBtn.disabled = false; clearBtn.disabled = false; generateBtn.addEventListener('click', processMoleculeAndGenerateKey); loadSmilesBtn.addEventListener('click', loadSmilesIntoComposer); clearBtn.addEventListener('click', clearWholeApp); } } function processMoleculeAndGenerateKey() { if (!rdkitModule || !composerApp) return; const chemObj = composerApp.getChemObj(); if (!chemObj || chemObj.isEmpty()) { inchikeyOutput.innerText = &quot;Error: Drawing board empty!&quot;; smilesOutput.value = &quot;&quot;; showStatus(&quot;Conversion failed: The drawing board is empty!&quot;, false); return; } try { const molfileData = Kekule.IO.saveFormatData(chemObj, 'mol'); if (!molfileData || molfileData.trim() === &quot;&quot;) { showStatus(&quot;Conversion failed: Structural generation error.&quot;, false); return; } const rdkitMol = rdkitModule.get_mol(molfileData); if (rdkitMol) { const canonicalSmiles = rdkitMol.get_smiles(); smilesOutput.value = canonicalSmiles; const inchi = rdkitMol.get_inchi(); if (inchi) { inchikeyOutput.innerText = rdkitModule.get_inchikey_for_inchi(inchi); showStatus(&quot;Success: 2D structure converted to identifiers!&quot;, true); } else { inchikeyOutput.innerText = &quot;Error: InChI failed.&quot;; showStatus(&quot;Partial conversion: InChI coordinates generation failed.&quot;, false); } rdkitMol.delete(); } else { inchikeyOutput.innerText = &quot;Error: Chemistry Validation Error.&quot;; showStatus(&quot;Conversion failed: Invalid chemical structure detected.&quot;, false); } } catch (error) { console.error(error); showStatus(&quot;Error: A critical parser exception occurred.&quot;, false); } } function loadSmilesIntoComposer() { if (!rdkitModule || !composerApp) return; const inputSmiles = smilesOutput.value.trim(); if (!inputSmiles || inputSmiles === &quot;Error: SMILES field is empty!&quot;) { showStatus(&quot;Conversion failed: Please enter a SMILES string.&quot;, false); return; } let rdkitMol = null; try { rdkitMol = rdkitModule.get_mol(inputSmiles); if (rdkitMol) { const molfile = rdkitMol.get_new_coords(); const chemObj = Kekule.IO.loadFormatData(molfile, 'mol'); if (chemObj) { composerApp.setChemObj(chemObj); const inchi = rdkitMol.get_inchi(); if (inchi) { inchiOutput.innerText = inchi; inchikeyOutput.innerText = rdkitModule.get_inchikey_for_inchi(inchi); } showStatus(&quot;Success: SMILES loaded to 2D canvas!&quot;, true); } else { showStatus(&quot;Conversion failed: Render engine error.&quot;, false); } rdkitMol.delete(); } else { showStatus(&quot;Conversion failed: Unparseable or invalid SMILES string.&quot;, false); } } catch (error) { console.error(error); showStatus(&quot;Error: A parser exception occurred handling the SMILES code.&quot;, false); if (rdkitMol) rdkitMol.delete(); } }</script></body></html>" title="2D Structure Tool"></iframe>
"""

# JS snippet for collapsible buttons
collapsible_js = """
<script>
    var coll = document.getElementsByClassName("collapsible");
    var i;
    for (i = 0; i < coll.length; i++) {
    coll[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var content = document.getElementsByClassName("inner")[0];
        var nextContent = this.nextElementSibling;
        console.log("search",content);
        console.log("sibling",nextContent);
        if (content.style.display === "block") {
        content.style.display = "none";
        } else {
        content.style.display = "block";
        }
    });
    }
    </script>   
"""

# Full collapsible + div for the editor
super_template = """{% block contents %}
    <div class="content">

    {{ super() }}
    
    <div class="divider"> </div>
    <button type="button" class="collapsible">Molecule editor to get InChIKeys for Locate Molecule (drop-down)</button>
    <div class="divider"> </div>
    <div class="inner">
    """ + iframe_editor + "</div> \n </div> \n" + collapsible_js + "{% endblock %}"

def cluster_nodes(descriptors, n_clusters="silhouettes", verbose=True):
    node_ids = list(descriptors.keys())
    X = np.array([descriptors[n] for n in node_ids])

    if n_clusters == "silhouettes":
        silhouettes = []
        k_values = range(2, int(np.sqrt(len(X) / 2)))  # candidate k values
        tolerance = 0.005
        window = 5
        for k in k_values:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10).fit(X)
            score = silhouette_score(X, kmeans.labels_)
            silhouettes.append(score)
            std = np.std(np.array(silhouettes))
            if verbose:
                print(f"Silhouette iteration {k} = {score:.3f}")
            if len(silhouettes) >= window:
                std = np.std(silhouettes[-window:])
                if std < tolerance:
                    print(f"Converged at k={k} (std={std:.4f})")
                    break
        # pick the k that gave the max silhouette
        best_index = np.argmax(silhouettes)
        n_clusters = k_values[best_index]  # <-- use k_values, not X
        if verbose:
            print("## Optimal number of clusters:", n_clusters)

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=n_clusters)
    cluster_labels = kmeans.fit_predict(X)

    clusters = {
        node_id: int(label) for node_id, label in zip(node_ids, cluster_labels)
    }
    return clusters


def assign_coordinates(graph, clusters):
    """
    Creates positions for NetworkX nodes such that nodes in the same cluster
    are closer together.
    """
    # place clusters on a circle, nodes in cluster randomly around center
    cluster_centers = {}
    n_clusters = len(set(clusters.values()))
    angle_step = 2 * np.pi / n_clusters

    for i, cluster_id in enumerate(sorted(set(clusters.values()))):
        cluster_centers[cluster_id] = np.array(
            [np.cos(i * angle_step), np.sin(i * angle_step)]
        )

    pos = {}
    for node in graph.nodes:
        center = cluster_centers[clusters[node]]
        # small random offset within cluster
        offset = np.random.normal(scale=0.1, size=2)
        pos[node] = center + offset
    return pos


def complete_cluster(node_list, cluster_dict, cluster_idx):
    """
    Function to assign a cluster index for unassigned nodes in a dictionary
    """
    unassigned = set(node_list).difference(set(cluster_dict.keys()))
    new_assignments = [(nd, cluster_idx) for nd in unassigned]
    cluster_dict.update(dict(new_assignments))
    return cluster_dict
    
def build_dashboard(G, compounds, title,outfile,size=(1400,800), 
                    layout_function="kamada_kawai",  
                    map_field="energy", verbose=True, add_editor=True,
                    **kwargs):
    """
    Wrapper function to generate HTML visualizations for a given network.

    Input:
    - G (nx.Graph): object as generated from RXReader. For profile support,
      it should contain a graph["pathList"] property.
    - title (str): title for the visualization.
    - outfile (str): name of the output HTML file.
    - size (tuple): tuple of integers, size of the visualization in pixels.
    - layout_function (nx.object, optional): Function to generate graph layout.
    - map_field (str): name of the field used for node coloring.

    Output:
    - lay (bokey.obj): Bokeh layout as generated by full_view_layout()
    """
    try:
        assert len(G.edges) > 0
    except AssertionError:
        print(
            "## WARNING! Chemical reaction network has no reactions. Check "
            "whether model variables (e.g., electronic method, solvent etc) "
            "fit model data in the MongoDB. Aborting html creation."
        )
        return None

    if verbose:
        print("## Writing {f1} output file".format(f1=outfile))
    # Define sizing
    w1 = int(size[0] * 4 / 7)
    w2 = int(size[0] * 3 / 7)
    wu = int(size[0] / 7)
    h = int(size[1] * 6 / 8)

    sizing_dict = {"w1": w1, "w2": w2, "wu": wu, "h": h}

    # Define custom classes
    jsmol_script = (
        "https://cdn.jsdelivr.net/gh/dgarayr/"
        "jsmol_to_bokeh/jsmol_to_bokeh.min.js"
    )
    style_template = """
    {% block postamble %}
	<script type="text/javascript" src="https://cdn.jsdelivr.net/gh/dgarayr/jsmol_to_bokeh/jsmol_to_bokeh.min.js"></script>
    <style>
    .bk-root .bk-btn-default {
        font-size: 1.2vh;
    }
    .bk-root .bk-input {
        font-size: 1.2vh;
        padding-bottom: 5px;
        padding-top: 5px;
    }
    .bk-root .bk {
        font-size: 1.2vh;
    }
    .bk-root .bk-clearfix{
        padding-bottom: 0.8vh;
    }
    .collapsible {
        background-color: #ffffff;
        cursor: pointer;
        padding: 18px;
        width: 100%;
        border: none;
        text-align: left;
        outline: none;
        font-size: 15px;
    }
    .active, .collapsible:hover {
        background-color: #cccccc;
        color: black;
        height: 50%;
    }
    .inner {
        padding: 0 18px;
        display: none;
        overflow: hidden;
        background-color: white;
    }
    .divider {
        height: 3px;
        width: 100%;
        background: #000000;
        cursor: col-resize;
    }
    </style>
    {% endblock %}
    """


    if add_editor:
        style_template += super_template 

    if layout_function == 'KMeans':  # custom clustering of nodes
        # this should be modifiable later
        cluster_property = "xyzdes"
        descriptors, prop_values = {}, []
        for nd in G.nodes(data=True):
            prop = nd[1][cluster_property]
            prop_values.append(prop)
        # prop_values,flag = aggregate_property(G,cluster_property)#,"none")
        # print("propvalues", prop_values)
        descriptors = dict(zip(G.nodes(), prop_values))
        clusters = cluster_nodes(descriptors)  # , n_clusters=5)
        posx = assign_coordinates(G, clusters)
    else:
        layout_function = getattr(nx, f"{layout_function}_layout")
        posx = layout_function(G)

    # Add model field to all nodes and edges & also vibrations
    arxviz.add_models(G)

    # Other properties
    node_size = kwargs.get("node_size", 30)

    # Bokeh-powered visualization via RXVisualizer
    bk_fig, bk_graph = arxviz.bokeh_network_view(
        G,
        positions=posx,
        graph_title=title,
        width=w1,
        height=h,
        map_field=map_field,
        hide_energy=True,
        palette=kwargs["palette"],
        qual_mapping=kwargs.get("qual_mapping", {}),
    )
    bk_graph.node_renderer.glyph.size = node_size

    # bk_graph.selection_policy = bkm.NodesAndLinkedEdges()
    bk_graph.selection_policy = bkm.EdgesAndLinkedNodes()

    # Modify the hovering tools to add additional fields
    # removing the previous ones first
    valid_tools = [tool for tool in bk_fig.tools if tool.description]
    old_hovers = [tool for tool in valid_tools if "hover" in tool.description]
    for tool in old_hovers:
        bk_fig.tools.remove(tool)

    # custom edge hovering to reduce noise
    #
    hover_edgeJS = """
    var erend = graph.edge_renderer.data_source
    var label1 = String.fromCharCode(916).concat("E1")
    var label2 = String.fromCharCode(916).concat("E2")
    if (cb_data.index.indices.length > 0) {
        var ndx = cb_data.index.indices[0]
        var tsname = erend.data["name"][ndx]
        if (tsname.includes('TSb')){
            hover.tooltips = [["tag","@name"]]
        } else {
            hover.tooltips = [["tag","@name"],["charge","@chargeStr"], ["pfcost", "@pfcostStr"],
                                ["multiplicity","@multiplicityStr"],["formula","@formulaStr"],
                                [label1,"@deltaE1"],[label2,"@deltaE2"]]
        }
    }
    """

    hide_barrlessJS = """
        var erend = graph.edge_renderer.data_source
        var nrend = graph.node_renderer.data_source
        var edgenames = erend.data["name"]
        var nodenames = nrend.data["name"]
        var numEdges = edgenames.length
        var numNodes = nodenames.length
        var statusCounter = counter[0]
        var connectedNodes = []
        var labsNodes = figure.center[2].source.data
        if (statusCounter == 0){
        // remove and set 1
            for (let j = 0; j < numEdges; j++){
                var is_tsb = edgenames[j].includes("TSb")
                if (is_tsb) {
                    erend.data["start"][j] = null
                    erend.data["end"][j] = null
                }
                else {
                    connectedNodes.push(erend.data["start"][j])
                    connectedNodes.push(erend.data["end"][j])
                }
            }
            for (let i = 0; i < numNodes; i++){
                var nname = nodenames[i]
                if (!connectedNodes.includes(nname)) {
                    nrend.data["index"][i] = null
                    labsNodes["nnames"][i] = " "
                }
            }
            statusCounter = 1
        } else {
        // restore and set counter back to zero
             for (let j = 0; j < numEdges; j++){
                var is_tsb = edgenames[j].includes("TSb")
                if (is_tsb){
                    erend.data["start"][j] = backupEdgeRoutes["start"][j]
                    erend.data["end"][j] = backupEdgeRoutes["end"][j]
                }
            }
            for (let i = 0; i < numNodes; i++) {
                if (nrend.data["index"][i] == null){
                    nrend.data["index"][i] = backupNodes["index"][i]
                    labsNodes["nnames"][i] = backupNodes["name"][i]
                }
            }
            statusCounter = 0
        }
        counter[0] = statusCounter
        nrend.change.emit()
        erend.change.emit()
        """

    # Custom locateMolecule function to support search by SMILES
    locateMolecule = """
        // source - source object for JSMol
        // pass graph and fetch node and edge renderers
        // from fig, we modify x_range and y_range.
        // Default plot starts from -1.2 to 1.2,
        var nrend = graph.node_renderer.data_source
        var erend = graph.edge_renderer.data_source
        var layout = graph.layout_provider.graph_layout
        // fetch the query in the data sources, choosing the appropiate
        // renderer depending on the query
        var mol_query = text_input.value
        if (mol_query.includes("TS") || mol_query.includes("ts")) {
            var renderer = erend
            var other_renderer = nrend
            var pool_smiles = []
            var pool_inchikeys = []
        } else {
            var renderer = nrend
            var other_renderer = erend
            var pool_smiles = renderer.data["smilesStr"]
            var pool_inchikeys = renderer.data["inchikey"]
        }
        var pool_names = renderer.data["name"]

        // split species joined by + sign, smiles by //,
        // inchikeys are already listd
        var pool_species = pool_names.map((name) => name.split("+"))
        var pool_smiles_split = pool_smiles.map((smi) => smi.split("//"))
        var split_ikey = function(ikey){
            var terms = [ikey.slice(0,14),ikey.slice(15,25),ikey];
            return terms}
        var pool_inchikeys_split = pool_inchikeys
                    .map((ikeys) => ikeys.map(split_ikey).flat())

        // function to match results in the array
        var getSubstringIndices = function(arr,query){
            return arr.reduce(
                    function(matches,tgt,i){
                            if (tgt.includes(query))
                        {matches.push(i)};
                            return matches;
                    },
                    []);
        }

        if (mol_query.includes("+")) {
            var ndx_u1 = pool_names.indexOf(mol_query)
            if (ndx_u1 < 0) {var ndx1 = []}
            else {var ndx1 = [ndx_u1]}
        } else {
            var ndx1 = getSubstringIndices(pool_species,mol_query)
        }

        var ndx2 = getSubstringIndices(pool_smiles_split,mol_query)
        var ndx3 = getSubstringIndices(pool_inchikeys_split,mol_query)
        // check all -> only choose the ones having matches => if several do,
        // order of preference is ikey/smiles/name

        if (ndx3.length > 0){
            var ndx = ndx3
        } else if (ndx2.length > 0) {
            var ndx = ndx2
        } else if (ndx1.length > 0) {
            var ndx = ndx1
        } else {
            var ndx = []
        }

        // locate positions of the node or of the nodes defining an edge
        if (mol_query.includes("TS") || mol_query.includes("ts")) {
            var n1 = renderer.data["start"][ndx]
            var n2 = renderer.data["end"][ndx]
            var pos1 = layout[n1]
            var pos2 = layout[n2]
            var positions = new Array(2)
            positions[0] = 0.5*(pos1[0]+pos2[0])
            positions[1] = 0.5*(pos1[1]+pos2[1])
        } else {
            var positions = layout[pool_names[ndx[0]]]
        }
        if (ndx.length > 0) {
            // clearing other sel. avoids problems for model loading sometimes
            other_renderer.selected.indices = []
            renderer.selected.indices = ndx
            fig.x_range.start = positions[0] - 0.5
            fig.x_range.end = positions[0] + 0.5
            fig.y_range.start = positions[1] - 0.5
            fig.y_range.end = positions[1] + 0.5
        }
        """

    # Callback for exporting the currently viewed nodes
    exportCurrent = """
        var nrend = graph.node_renderer.data_source
        var nodeIndices = nrend.data["index"]
        var presentNodes = nodeIndices.filter((node) => node != null)
        var outDict = {nodes: presentNodes}
        // downloading a file
        function download(content, fileName, contentType) {
            var a = document.createElement("a")
            var file = new Blob([content], {type: contentType})
            a.href = URL.createObjectURL(file)
            a.download = fileName
            a.click()
        }
        download(
                    JSON.stringify(outDict),
                    "vizchemoton_export_sel.json",
                    "text/plain"
                )
    """

    tooltips = [
        ("tag", "@name"),
        ("charge", "@chargeStr"),
        ("multiplicity", "@multiplicityStr"),
        ("formula", "@formulaStr"),
        ("smiles", "@smilesStr"),
        ("pfcost", "@pfcostStr"),
    ]
    tooltips += kwargs.get("custom_hovers", [])

    hover_node = bkm.HoverTool(
        description="Node hover",
        renderers=[bk_graph.node_renderer],
        tooltips=tooltips,
        formatters={"@energy": "printf"},
    )
    bk_fig.add_tools(hover_node)
    hover_edge = bkm.HoverTool(
        description="Edge hover",
        renderers=[bk_graph.edge_renderer],
        formatters={"@energy": "printf"},
        line_policy="interp",
    )
    hover_edge.callback = bkm.CustomJS(
        args={"hover": hover_edge, "graph": bk_graph}, code=hover_edgeJS
    )
    bk_fig.add_tools(hover_edge)

    highl_callback = bkm.CustomJS(
        args={"graph": bk_graph},
        code=arxviz.js_callback_dict["highlightNeighbors"],
    )

    # We need edge backups
    edgesource = bk_graph.edge_renderer.data_source
    backup_edges = {
        "start": copy.deepcopy(edgesource.data["start"]),
        "end": copy.deepcopy(edgesource.data["end"]),
    }
    backup_nodes = {
        "index": bk_graph.node_renderer.data_source.data["index"],
        "name": bk_graph.node_renderer.data_source.data["name"],
    }
    hide_barrless_callback = bkm.CustomJS(
        args={
            "graph": bk_graph,
            "figure": bk_fig,
            "counter": [0],
            "backupNodes": backup_nodes,
            "backupEdgeRoutes": backup_edges,
        },
        code=hide_barrlessJS,
    )

    alt_ref_e = kwargs.get("alt_ref_energy",0.0)
    lay = arxviz.full_view_layout(bk_fig, bk_graph, G=G, sizing_dict=sizing_dict,alt_ref_energy=alt_ref_e)

    # add a button to the layout
    b_highlight = bkm.Button(
        label="Highlight neighbors", max_width=int(w1 / 6), align="center"
    )
    b_highlight.js_on_click(highl_callback)
    b_hidebarrless = bkm.Button(
        label="Hide barrierless", max_width=int(w1 / 6), align="center"
    )
    b_hidebarrless.js_on_click(hide_barrless_callback)

    sel_row = lay.children[0][0].children[2]
    sel_row.children[1].max_width = int(w1 / 6)
    sel_row.children = (
        sel_row.children[0:2]
        + [b_highlight, b_hidebarrless]
        + [sel_row.children[-1]]
    )

    # Export functionality
    export_callback = bkm.CustomJS(
        args={"figure": bk_fig, "graph": bk_graph}, code=exportCurrent
    )
    # Additional upper button -> readjust spacing to fit
    up_row = lay.children[0][0].children[0]
    b_export = bkm.Button(
        label="Export current nodes", max_width=int(w1 / 6), align="center"
    )
    b_export.js_on_click(export_callback)

    for item in up_row.children:
        item.max_width = int(w1 / 6)

    up_row.children.append(b_export)

    # Modify the callback of the locate molecule button
    text_input = sel_row.children[0]
    js_mol_locator_nw = bkm.CustomJS(
        args={"graph": bk_graph, "fig": bk_fig, "text_input": text_input},
        code=locateMolecule,
    )
    sel_button = sel_row.children[1]
    sel_button.js_event_callbacks["button_click"] = [js_mol_locator_nw]
    sel_button.js_on_click(js_mol_locator_nw)

    bokeh.plotting.output_file(outfile, title=title, mode="cdn")
    bokeh.plotting.save(lay, template=style_template)

    return lay, bk_fig, bk_graph


def scale_xyz_list(xyz, displ_vector=np.zeros(3)):
    """
    Bohr-to-angstrom scaling of a list of XYZ coordinates of the form
    [atom, [x, y, z]].

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with
    atom being a string and x,y,z floats.
    - displ_vector (np.Array, optional): for translating the geometry.

    Output:
    - xyz_nw (list): scaled XYZ coordinates in the same format as the input.
    """

    bohr_to_ang = 0.529
    xyz_arr = np.array([item[1] for item in xyz]) * bohr_to_ang + displ_vector
    xyz_nw = [[item[0], list(xyz_arr[ii])] for ii, item in enumerate(xyz)]
    return xyz_nw


def xyz_list_to_xyz_block(xyz):
    """
    Transform a list of xyz coordinates [atom, [x, y, z]] into a string block.

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with
    atom being a string and x,y,z floats.

    Output:
    - xyz_block (str): newline-joined block of the form
    a1,x1,y1,z2\na2,x2,y2,z2...
    """

    xyz_block = "\n".join(
        ["%s %.6f %.6f %.6f" % (item[0], *item[1]) for item in xyz]
    )
    return xyz_block


def formula_from_xyz_block(xyz):
    """
    Generates the molecular formula for a given XYZ geometry.

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with
    atom being a string and x,y,z floats.

    Output:
    - formula (str): molecular formula from the input geometry.
    """
    labels = [item[0] for item in xyz]
    counter_list = sorted(Counter(labels).items())
    formula = ""
    for atom, ct in counter_list:
        if ct == 1:
            formula += atom
        else:
            formula += "%s%d" % (atom, ct)
    return formula


def sort_edge_names(edge_tuple):
    """
    Helper function to sort edge tuples lexicographically.

    Input:
    - edge_tuple (tuple): edge specification as a pair of node names.

    Output:
    - lexico_tuple (tuple): lexicographically sorted tuple.
    """
    n1, n2 = [int(nd) for nd in edge_tuple]
    srt_pair = sorted([n1, n2])
    lexico_tuple = tuple((str(nd) for nd in srt_pair))
    return lexico_tuple


def preprocess_compounds(compounds):
    """
    Helper function to process compounds properties.
    """
    tgt_vars = ["energy", "charge", "multiplicity", "pfcost"]
    for comp in compounds.values():
        for vv in tgt_vars:
            if not isinstance(comp[vv], list):
                comp[vv] = [comp[vv]]


def build_graph_edges(reaction_list):
    """
    Build graph edges.
    """
    return [(it[0], it[1], {"tsidx": it[2]}) for it in reaction_list]


def get_node_name_and_geometry(comp, dist_adduct, bohr_to_ang):
    """
    Helper function to retrieve the node name and geometry.
    """
    if isinstance(comp["crn_id"], list):
        comp["crn_id"] = "+".join(comp["crn_id"])
    if "+" in comp["crn_id"]:  # or isinstance(comp["crn_id"], list):
        # node_name = "+".join(comp["crn_id"])
        node_name = comp["crn_id"]
        xyz_list = comp["xyz"]
        xyz0_arr = np.array([item[1] for item in xyz_list[0]]) * bohr_to_ang
        cntr = xyz0_arr.mean(axis=0)
        xyz0 = [
            [item[0], list(xyz0_arr[ii])]
            for ii, item in enumerate(xyz_list[0])
        ]
        xyz_full = xyz0

        for ii, xyz in enumerate(xyz_list[1:]):
            displ_vec = cntr + (ii + 1) * dist_adduct
            xyz_arr = np.array([it[1] for it in xyz]) * bohr_to_ang + displ_vec
            xyz_nw = [
                [item[0], list(xyz_arr[ii])] for ii, item in enumerate(xyz)
            ]
            xyz_full += xyz_nw
    else:
        node_name = comp["crn_id"]
        xyz_list = [comp["xyz"]]
        xyz_arr = np.array([item[1] for item in xyz_list[0]]) * bohr_to_ang
        xyz_full = [
            [item[0], list(xyz_arr[ii])] for ii, item in enumerate(xyz_list[0])
        ]

    return node_name, xyz_full, xyz_list


def add_node_attributes(
    graph, compounds, node_renaming, dist_adduct, bohr_to_ang
):
    """
    Add nodes attributes to the graph for building the HTML file.
    """
    for nd in graph.nodes(data=True):
        comp = compounds[nd[0]]
        tmp = get_node_name_and_geometry(comp, dist_adduct, bohr_to_ang)
        node_name, xyz_full, xyz_list = tmp
        node_renaming[nd[0]] = node_name
        strtmp = "%s %.6f %.6f %.6f"
        _xyz = "\n".join([strtmp % (item[0], *item[1]) for item in xyz_full])
        nd[1]["geometry"] = _xyz
        nd[1]["energy"] = sum(comp["energy"])
        nd[1]["ZPVE"] = 0.0
        nd[1]["name"] = node_name
        nd[1]["degree"] = graph.degree(nd[0])
        nd[1]["charge"] = comp["charge"]
        nd[1]["multiplicity"] = comp["multiplicity"]
        nd[1]["formula"] = [formula_from_xyz_block(xyz) for xyz in xyz_list]
        nd[1]["neighbors"] = list(graph.neighbors(nd[0]))
        nd[1]["smiles"] = str(comp.get("smiles", "None")).split("//")
        nd[1]["inchikey"] = str(comp.get("inchikey", "None")).split("//")
        nd[1]["xyzdes"] = comp["xyzdes"]
        nd[1]["pfcost"] = comp["pfcost"]


def add_edge_attributes(graph, compounds):
    """
    Add edges attributes to the graph for building the HTML file.
    """
    compounds_renamed = {}
    for _, comp in compounds.items():
        if isinstance(comp["crn_id"], list):
            new_key = "+".join(comp["crn_id"])
            compounds_renamed[new_key] = comp
        else:
            compounds_renamed[comp["crn_id"]] = comp
    for ii, ed in enumerate(graph.edges(data=True)):
        e1, e2 = [sum(compounds_renamed[nd]["energy"]) for nd in ed[0:2]]
        if ed[2]["tsidx"] is None or ed[2]["tsidx"] == "None":
            e_ts = max(e1, e2)
            ed[2]["name"] = "TSb_%04d" % ii
            ed[2]["geometry"] = None
            ed[2]["energy"] = 0.0
            ed[2]["ZPVE"] = 0.0
            delta_e1 = (e_ts - e1, ed[0])
            delta_e2 = (e_ts - e2, ed[1])
            ed[2]["deltaE1"] = "%.2f (%s)" % delta_e1
            ed[2]["deltaE2"] = "%.2f (%s)" % delta_e2
            continue
        ts_compound = compounds[ed[2]["tsidx"]]
        xyz_list = [ts_compound["xyz"]]
        geom = scale_xyz_list(xyz_list[0])
        ed[2]["geometry"] = xyz_list_to_xyz_block(geom)
        ed[2]["name"] = ts_compound["crn_id"]
        e_ts = sum(ts_compound["energy"])
        delta_e1 = (e_ts - e1, ed[0])
        delta_e2 = (e_ts - e2, ed[1])
        ed[2]["deltaE1"] = "%.2f (%s)" % delta_e1
        ed[2]["deltaE2"] = "%.2f (%s)" % delta_e2
        ed[2]["energy"] = e_ts
        ed[2]["ZPVE"] = 0.0
        ed[2]["charge"] = ts_compound["charge"]
        ed[2]["multiplicity"] = ts_compound["multiplicity"]
        ed[2]["formula"] = [formula_from_xyz_block(xyz) for xyz in xyz_list]


def format_string_attributes(graph):
    """
    Processes node & edge attributes that are shown as strings in the
    final dashboard
    """
    node_attrs = ["charge", "multiplicity", "formula", "smiles", "pfcost"]
    edge_attrs = ["charge", "multiplicity", "formula"]
    for nd in graph.nodes(data=True):
        for tgt in node_attrs:
            nd[1][tgt + "Str"] = "//".join([str(item) for item in nd[1][tgt]])

    for ed in graph.edges(data=True):
        for tgt in edge_attrs:
            if tgt not in ed[2].keys():
                continue
            ed[2][tgt + "Str"] = "//".join([str(item) for item in ed[2][tgt]])

    return None


def _clean_empty_entries(compounds):
    """
    Preprocessing function to remove empty dictionaries in compounds object.
    """
    filtered_compounds = {}
    for c in compounds:
        if compounds[c] != {}:
            filtered_compounds[c] = compounds[c]
    return filtered_compounds


def process_graph(reaction_list, compounds, dist_adduct=3.0):
    """
    Wrapper function to generate a nx.Graph from a list of reactions and a
    dictionary of compounds, including XYZ-formatted geometries where
    individual geometries of the species forming adducts are joined.
    """
    compounds = _clean_empty_entries(compounds)
    bohr_to_ang = 0.529177
    graph = nx.Graph()
    edge_list = build_graph_edges(reaction_list)
    graph.add_edges_from(edge_list)
    node_renaming = {}
    preprocess_compounds(compounds)
    add_node_attributes(
        graph, compounds, node_renaming, dist_adduct, bohr_to_ang
    )
    nx.relabel_nodes(graph, node_renaming, copy=False)

    # update neighbors after renaming
    for nd in graph.nodes(data=True):
        nd[1]["neighbors"] = list(graph.neighbors(nd[0]))

    add_edge_attributes(graph, compounds)
    format_string_attributes(graph)

    return graph

def property_list_flatter(prop_dict,sep="//"):
    """For a list of node/edge properties, collapse lists into strings (if no propertyStr property exists yet)
    (and stringify None)"""
    to_add = {}
    to_remove = []
    for k,v in prop_dict.items():
        if v is None:
            prop_dict[k] = "None"
            continue 
        if not(isinstance(v,list)):
            continue 
        if f"{k}Str" not in prop_dict.keys():
            to_add[f"{k}Str"] = sep.join([str(val) for val in v])
        to_remove.append(k)
    for k in to_remove:   
        del prop_dict[k]
    prop_dict.update(to_add)
    return prop_dict
            
def save_graph(graph,filename):
    """
    Wrapper function to save the graph to GraphML format (Gephi-compatible). Must collapse lists into strings.
    """
    gwork = graph.copy()
    for nd in gwork.nodes(data=True):
        property_list_flatter(nd[1]) 
    for ed in gwork.edges(data=True):
        property_list_flatter(ed[2])
    nx.write_graphml(gwork,path=filename)
    return None

def format_value_list(val_list, fmt="%.4f", sep="//"):
    """
    TO-DO
    """
    fvlist = [fmt % vv if vv is not None else "None" for vv in val_list]
    return sep.join(fvlist)


def aggregate_property(Gx, prop_name, agg_func="mean", na_value=0):
    """
    TO-DO
    """
    fmap = {
        "max": np.max,
        "min": np.min,
        "mean": np.mean,
        "sum": np.sum,
        "none": lambda x: x,
    }
    func = fmap.get(agg_func, np.mean)
    agg_values = []
    flags = []
    for nd in Gx.nodes(data=True):
        prop = nd[1][prop_name]
        flag = 0
        if isinstance(prop, float) or isinstance(prop, int):
            val = prop
        elif isinstance(prop, list):
            values = [
                item if item is not None else np.nan
                for item in nd[1][prop_name]
            ]
            mask = np.isnan(values)
            values = np.where(mask, na_value, values)
            if np.all(mask):
                flag = 2
            elif np.any(mask):
                flag = 1
            val = func(values)
        elif prop is None:
            val = na_value
            flag = 2

        agg_values.append(val)
        flags.append(flag)
    return agg_values, flags


### Path management functions - July 2026
def identify_balanced_node(graph,node,neighbor):
    """Convenience function to locate a suitable balanced (A+B) node in a graph, 
    given one of the involved compounds (A or B) and a neighboring node"""
    if "ts" in neighbor.lower():
        # Locate by name
        pair = [ed[0:2] for ed in graph.edges(data="name") if ed[2] == neighbor][0]
        onode = [item for item in pair if item != neighbor][0]
    else:
        # Get all neighbors of the neighbor and filter by name
        sel = [nd for nd in graph[neighbor] if node in nd]
        onode = sel[0]
    return onode

def path_adjuster(graph,path_list_raw):
    """Processes existing paths (e.g. from pathfinder) to make them compliant with amk-tools considerations: balanced nodes 
    (for start and end) and uppercase TS labels"""
    out_path_list = []
    for path in path_list_raw:
        fmt_path = [entry.replace("ts","TS") if "ts" in entry else entry for entry in path[1:-1]]
        # processing first and last entries
        source = path[0]
        target = path[-1]
    
        # if they are a node, proceed
        if source in graph.nodes():
            fmt_path = [source] + fmt_path
        else:
            neigh = path[1]
            onode = identify_balanced_node(graph,source,neigh)
            fmt_path = [onode] + fmt_path
            
        if target in graph.nodes():
            fmt_path += [target]
        else:
            neigh = path[-2]
            onode = identify_balanced_node(graph,target,neigh)
            fmt_path += [onode]

        out_path_list.append(fmt_path)
    return out_path_list

def prepare_extrema(graph,node):
    """For a given string defining an entity in the CRN, retrieves the corresponding node (if existing)
    or the collection of all the nodes where that entity participates (e.g. A -> A+B, A+C, D+A...)"""
    if node not in graph.nodes():
        nodeset = [nd for nd in graph.nodes if node in nd.split("+")]
    else:
        nodeset = [node]
    return nodeset

def check_pfcosts_in_path(graph,path,limit_value=9999):
    """For a given path composed of valid nodes through the graph, returns a list with all
    the node costs through the path."""
    path_nodes = [item for item in path if "ts" not in item.lower()]
    costs = [sum(graph.nodes[nd].get('pfcost',[limit_value])) for nd in path_nodes]
    return costs

def generate_paths(graph,source,target,max_length=6,check_costs=False,Npaths_filt=10):
    """Wrapper function to generate valid paths for a given graph, using a maximum cutoff length (longer paths
    will not be considered to limit path search cost). Source and target are automatically checked for node 
    validity, finding combinations through prepare_extrema. The pfcost property of nodes can be used to 
    accumulate the cost of each path and rank them"""
    src_nodes = prepare_extrema(graph,source)
    end_nodes = prepare_extrema(graph,target)
    found_paths = arx.add_paths(graph,src_nodes,end_nodes,cutoff=max_length)
    found_paths_cln = path_adjuster(graph,found_paths)

    if check_costs:
        all_costs = [(idx,sum(check_pfcosts_in_path(graph,path))) for idx,path in enumerate(found_paths_cln)]
        srt_costs = sorted(all_costs,key=itemgetter(1))
        Np = min(len(found_paths_cln),Npaths_filt)
        sel_paths = [found_paths_cln[idx] for idx,cost in srt_costs[0:Np]] 
    else:
        sel_paths = found_paths_cln 

    return sel_paths

def get_energy_ref(graph,path_list):
    """Gets the energy of all starting nodes involved in the path list and returns the lowest-energy one to
    be used as reference for relative energies"""
    start_nodes = [path[0] for path in path_list]
    all_energies = [graph.nodes[nd]["energy"] for nd in start_nodes]
    idx = np.argmin(all_energies)
    return (start_nodes[idx],all_energies[idx])