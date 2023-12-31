import exporter from "./kekule.moduleEnvInits.esm.mjs";
import "./kekule.js";
import "./mins/root.min.js";
import "./mins/localization.min.js";
import "./mins/localizationData.zh.min.js";
import "./mins/common.min.js";
import "./mins/core.min.js";
import "./mins/html.min.js";
import "./mins/io.min.js";
import "./mins/render.min.js";
import "./mins/spectroscopy.min.js";
import "./mins/widget.min.js";
import "./mins/chemWidget.min.js";
import "./mins/webComponent.min.js";
import "./mins/algorithm.min.js";
import "./mins/calculation.min.js";
import "./mins/data.min.js";
import "./mins/emscripten.min.js";
import "./mins/openbabel.min.js";
import "./mins/indigo.min.js";
import "./mins/inchi.min.js";
let { Kekule, Class, ClassEx, ObjectEx, DataType} = exporter();
export { Kekule, Class, ClassEx, ObjectEx, DataType};
if(!Kekule.scriptSrcInfo.modules)Kekule.scriptSrcInfo.modules=[];
Kekule.ArrayUtils.pushUnique(Kekule.scriptSrcInfo.modules, ["lan", "root", "localization", "localizationData", "localizationData.zh", "common", "core", "html", "io", "render", "spectroscopy", "widget", "chemWidget", "webComponent", "algorithm", "calculation", "data", "emscripten", "openbabel", "indigo", "inchi"]);
if (typeof(Kekule) !== 'undefined') { Kekule._loaded(); }