/*
Copyright (c) 2016 - 2024, Adrian Dusa
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, in whole or in part, are permitted provided that the
following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The names of its contributors may NOT be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ADRIAN DUSA BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

$( function() {
$('#main_menu').smartmenus({
    subMenusSubOffsetX: 6,
    subMenusSubOffsetY: -8
});
var spinnerstarted = false;
function positionSpinner() {
    if ($("#background").length) {
        var result = document.getElementById("result_main");
        $("#background").css({
            "width": $("#result").width(),
            "height": $("#result").height() - 20
        });
        $("#background").position({my: "left top", at: "left top", of: result, collision: "none"});
        $("#spinner").position({my: "left top", at: "left+" + ($("#result").width() / 2 - 40) + "px top+" + ((($("#result").height() - 20) / 2) - 40) + "px", of: result, collision: "none"});
    }
}
function createSpinner() {
    var background = document.createElement("div");
    background.id = "background";
    document.body.appendChild(background);
    $("#background").addClass("spinnerbackground");
    background = document.getElementById("background");
    var spinner = document.createElement("div");
    spinner.id = "spinner";
    background.appendChild(spinner);
    $("#spinner").addClass("cssload-loader");
    spinner = document.getElementById("spinner");
    var one = document.createElement("div");
    one.id = "one";
    spinner.appendChild(one);
    $("#one").addClass("cssload-inner cssload-one");
    var two = document.createElement("div");
    two.id = "two";
    spinner.appendChild(two);
    $("#two").addClass("cssload-inner cssload-two");
    var three = document.createElement("div");
    three.id = "three";
    spinner.appendChild(three);
    $("#three").addClass("cssload-inner cssload-three");
    $("#background").css("z-index", "9099");
    spinnerstarted = true;
    positionSpinner();
}
function removeSpinner() {
    $("#background").remove();
    spinnerstarted = false;
}
var responseR = false;
var divid = 0;
var history = ["library(QCA)"];
var histindex = 1;
var activecommand = "";
var help = 0;
var visiblerows = 16, visiblecols = 7;
var gridset, datacover;
var deheight = 400, dewidth = 659;
var tempdirfile = "";
var dirfile = "";
var dirfile_chosen = ["dir", "~", ""];
var rects_width;
var canvas_height;
var tempdatainfo = {ncols: 0, nrows: 0, colnames: [], rownames: []};
var info = {data: null, tt: null, qmc: null};
var ovBox, input; 
var tastaRcommand = "";
var tasta = "enter";
var outres;
function reset_outres() {
    outres = {
        error: null,
        warning: null,
        output: null,
        message: null,
        tt: null,
        infobjs: null
    }
}
reset_outres();
var Rcommand = {
    "counter": 0,
    "command": "",
    "complete": "",
    "brackets": [],
    "thinfo": {},
    "scrollvh": {},
    "calibrate": {},
    "xyplot": {},
    "plotsize": [550/96, 550/96]
}
var scrollobj = {
    "counter": 0,
    "scrollvh": [],
    "dataset": ""
}
var scrolleftop = {};
var closeplot = 0;
var plotopen = false;
var tempcommand = "";
var txtcommand = "";
var escapeopen = false;
var changes = 0; 
var windowHeight = window.innerHeight;
var commandHeight = 100, resultHeight = 600;
if (windowHeight < 600) {
    resultHeight = 300;
    commandHeight = 75;
}
var string_command = "";
var dirfilevisit = false;
var minimize2R, tt2R, calib2R;
var colclicks = new Object;
var current_command = "";
var objname = "";
var testX = 80, testY = 33;
var scrollbarsWH = getScrollBarWidth();
var visibledata, updatecounter = 0, updatecounter2 = 0;
var pathcopy;
var outputcopy, coordscopy;
var dirsfilescopy;
var dirfilist = {
    refresh: true,
    value: 0
}
var thinfo = {
    "counter": 0,
    "dataset": "",
    "condition": "",
    "findth": false,
    "nth": 1
};
var ths = new Array(6);
var thsetter_content; 
var poinths;
function thinfoToR() {
    updatecounter = 0;
    responseR = false;
    thinfo.counter = 1 - thinfo.counter; 
    thinfo.dataset = commobj["calibrate"]["dataset"];
    thinfo.x = commobj["calibrate"]["x"];
    thinfo.findth = commobj["calibrate"]["findth"];
    thinfo.nth = commobj["calibrate"]["nth"];
    Shiny.onInputChange("thinfo", thinfo);
    doWhenDataPointsAreReturned();
}
function resetPoinths() {
    poinths = {
        "message": "",
        "vals": new Array(),
        "thvals": new Array(),
        "fuzzyvals": new Array(),
        "prettyx": new Array()
    }
}
resetPoinths();
var handles = new Array(6);
var lastvals = new Array();
var xyplotdata = new Array();
var rloadcycles = 0;
var ttfromR = {};
var papers = {}; 
var packages = {}; 
var cpl = { 
    command: 76,
    result: 74
}
var focus = "result";
var commobj = {};
var settings = {
    load: {
        name:      "load",
        title:     "Load data from attached packages",
        position:  {my: "left top", at: "left+5px top+33px", of: window, collision: "none"},
        resizable: false,
        width:     420,
        height:    275,
        inside: {
            package: {border: true, left: 15, top: 50, width: 160, height: 100},
            pkgdata: {border: true, left: 185, top: 50, width: 220, height: 100, dependency: "package"},
            title:   {border: true, left: 15, top: 160, width: 390, height: 60}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "package": "",
                "pkgdata": ""
            }
        }
    },
    import: {
        name:      "import",
        title:     "Import from text file",
        position:  {my: "left top", at: "left+5px top+33px", of: window, collision: "none"},
        resizable: false,
        width:     645,
        height:    433,
        inside:    {
            cols: {border: true, left:  14, top: 260, width: 225, height: 90},
            path: {border: true, left: 254, top:  62, width: 375, height:  40},
            dirs: {border: true, left: 254, top:  80, width: 375, height: 300}
        },
        reset: function(x) {
            x = (x == "import")?"read_table":x;
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "stdir": "",
                "objname": "",
                "sep": ",",
                "dec": ".",
                "header": true,
                "row_names": "",
                "nameit": true,
                "customname": false,
                "fileEncoding": "UTF-8-BOM"
            }
        }
    },
    export: {
        name:      "export",
        title:     "Export to text file",
        position:  {my: "left top", at: "left+5px top+33px", of: window, collision: "none"},
        resizable: false,
        width:     600,
        height:    380,
        inside: {
            dataset: {border: true, left:  18, top: 243, width: 195, height:  80},
            path:    {border: true, left: 230, top:  65, width: 355, height:  40},
            dirs:    {border: true, left: 230, top:  83, width: 355, height: 240}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "dataset": "", 
                "filename": "",
                "sep": ",",
                "dec": ".",
                "header": true,
                "caseid": "cases",
                "newfile": false
            };
        }
    },
    data_editor: {
        name:      "data_editor",
        title:     "Data editor",
        position:  {my: "left top", at: "left+" + testX + "px top+" + testY + "px", of: window, collision: "none"},
        resizable: true,
        width:     200, 
        height:    150, 
        inside:    {
            topleft:  {border: false, left:  0, top: 20, width: 70, height: 20},
            colnames: {border: false, left: 70, top: 20, width: (visiblecols + 1)*70, height: 20},
            rownames: {border: false, left:  0, top: 40, width: 70, height: (visiblerows + 1)*20},
            body:     {border: false, left: 70, top: 40, width: (visiblecols + 1)*70, height: (visiblerows + 1)*20}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = new Object;
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "dataset": "" 
            }
        }
    },
    calibrate: {
        name:      "calibrate",
        title:     "Calibrate",
        position:  {my: "left top", at: "left+80px top+33px", of: window, collision: "none"},
        resizable: true,
        width:     500,
        height:    390,
        minWidth:  500,
        minHeight: 390,
        inside:    {
            dataset: {border: true, left: 15, top: 52, width: 120, height: 100},
            x:       {border: true, left: 15, top: 182, width: 120, height: 140, dependency: "dataset"}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "dataset": "", 
                "x": "",
                "type": "fuzzy",
                "thresholds": ["", "", "", "", "", ""],
                "thnames": ["", "", "", "", "", ""],
                "thscopycrp": ["", "", "", "", "", ""],
                "thscopyfuzs": ["", "", "", "", "", ""], 
                "thscopyfuzb": ["", "", "", "", "", ""], 
                "logistic": true,
                "idm": "0.95",
                "ecdf": false,
                "below": "1",
                "above": "1",
                "same": true,
                "newvar": "",
                "increasing": true,
                "end": true,
                "findth": false,
                "nth": 1,
                "thsetter": false,
                "scrollvh": {},
                "jitter": false
            };
        }
    },
    recode: {
        name:      "recode",
        title:     "Recode",
        position:  {my: "left top", at: "left+80px top+33px", of: window, collision: "none"},
        resizable: false,
        width:     515,
        height:    340,
        inside: {
            dataset: {border: true, left: 14, top:  52, width: 120, height: 80},
            x:       {border: true, left: 14, top:  162, width: 120, height: 120, dependency: "dataset"},
            rules:   {border: true, left: 325, top: 202, width: 170, height: 80}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "dataset": "", 
                "x": "",
                "same": true,
                "newvar": "",
                "oldv": new Array(),
                "newv": new Array(),
                "scrollvh": {}
            };
        }
    },
    findRows: {
        name:      "findRows",
        title:     "Find incoherent configurations",
        position:  {my: "left top", at: "left+170px top+33px", of: window, collision: "none"},
        resizable: false,
        width:     464,
        height:    222,
        inside: {
            tt: {border: true, left:  14, top:  47, width: 138, height: 120}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "objname": "", 
                "nameit": false,
                "obj": "",
                "expression": "",
                "remainders": true,
                "type": [false, true, false, false],
                "remainders": true
            };
        }
    },
    tt: {
        name:      "tt",
        title:     "Truth table",
        position:  {my: "left top", at: "left+170px top+33px", of: window, collision: "none"},
        resizable: false,
        width:     464,
        height:    400,
        inside: {
            dataset:    {border: true, left:  14, top:  47, width: 138, height: 120},
            outcome:    {border: true, left: 162, top:  47, width: 138, height: 120, dependency: "dataset"},
            conditions: {border: true, left: 310, top:  47, width: 138, height: 120, dependency: "dataset"}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "objname": "", 
                "nameit": false,
                "dataset": "", 
                "outcome": new Array(1),
                "conditions": new Array(),
                "neg_out": false,
                "ic1": "1",
                "ic0": "",
                "pri_cut": "",
                "n_cut": "1",
                "exclude": "",
                "complete": false,
                "show_cases": false,
                "dcc": false,
                "sort_by": {"out": true, "incl": true, "n": true},
                "sort_sel": {"out": false, "incl": false, "n": false},
                "decreasing": true,
                "use_letters": false,
                "inf_test": "",
                "PRI": true
            };
        }
    },
    minimize: {
        name:       "minimize",
        title:      "Quine-McCluskey minimization",
        position:   {my: "left top", at: "left+170px top+33px", of: window, collision: "none"},
        resizable:  false,
        width:      464,
        height:     470,
        inside: {
            dataset:    {border: true, left:  14, top:  47, width: 138, height: 120},
            outcome:    {border: true, left: 162, top:  47, width: 138, height: 120, dependency: "dataset"},
            conditions: {border: true, left: 310, top:  47, width: 138, height: 120, dependency: "dataset"}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "objname": "", 
                "nameit": false,
                "dataset": "", 
                "outcome": new Array(),
                "conditions": new Array(),
                "neg_out": false,
                "ic1": "1",
                "ic0": "",
                "pri_cut": "",
                "n_cut": "1",
                "include": ["", ""],
                "exclude": "",
                "details": false,
                "row_dom": false,
                "all_sol": false,
                "first_min": false,
                "dir_exp": "",
                "inf_test": "",
                "pi_cons": "0",
                "pi_depth": "0",
                "sol_cons": "0",
                "sol_cov": "1",
                "sol_depth": "0",
                "use_letters": false,
                "PRI": true,
                "source": "tt" 
            };
        }
    },
    xyplot: {
        name:      "xyplot",
        title:     "XY plot",
        position:  {my: "left top", at: "left+240px top+33px", of: window, collision: "none"},
        resizable: true,
        width:     720,
        height:    567,
        minWidth:  720,
        minHeight: 567,
        inside: {
            dataset: {border: true, left: 13, top:   50, width: 150, height: 80},
            x:       {border: true, left: 13, top:  165, width: 150, height: 80, dependency: "dataset"},
            y:       {border: true, left: 13, top:  280, width: 150, height: 80, dependency: "dataset"}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "dataset": "", 
                "x": "",
                "y": "",
                "sufnec": "sufficiency",
                "pof": true,
                "mdguides": true,
                "labels": false,
                "fill": true,
                "jitter": false,
                "negy": false,
                "negx": false
            };
            xyplotdata = new Array();
        }
    },
    venn: {
        name:      "venn",
        title:     "Venn diagram",
        position:  {my: "left top", at: "left+240px top+33px", of: window, collision: "none"},
        resizable: true,
        width:     430,
        height:    500,
        minWidth:  430,
        minHeight: 500
    },
    saveRplot: {
        name:       "saveRplot",
        title:      "Save R plot window",
        position:   {my: "left top", at: "left+5px top+33px", of: window, collision: "none"},
        resizable:  false,
        width:      630,
        height:     400,
        inside: {
            preview: {border: false, left:   0, top:  20, width: 240, height: 250},
            path:    {border: true, left: 240, top:  62, width: 375, height:  40},
            dirs:    {border: true,  left: 240, top:  80, width: 375, height: 260}
        },
        reset: function(x) {
            if (commobj[x] === void 0) {
                commobj[x] = {};
            }
            commobj[x] = {
                "counter": (commobj[x].counter === void 0)?0:commobj[x].counter,
                "filename": "",
                "type": "png"
            };
        }
    },
    about: {
        name:       "about",
        title:      "About this software",
        position:   {my: "left top", at: "left+300px top+33px", of: window, collision: "none"},
        resizable:  false,
        width:      455,
        height:     440 + 10*(navigator.browserType == "Firefox")
    },
    changes: {
        name:       "changes",
        title:      "Change log",
        position:   {my: "left top", at: "left+300px top+33px", of: window, collision: "none"},
        resizable:  false,
        width:      550,
        height:     420 + 10*(navigator.browserType == "Firefox")
    },
    command: {
        name:      "command",
        title:     "Command constructor",
        position:  {my: "right top", at: "right-10px top+33px", of: window, collision: "fitflip"},
        resizable: true,
        width:     665, 
        height:    commandHeight, 
        closable: false
    },
    result: {
        name:      "result",
        title:     "R console",
        position:  {my: "left top", at: "left bottom+4px", of: "#command", collision: "fitflip"},
        resizable: true,
        width:     665, 
        height:    resultHeight, 
        closable: false
    },
    plotdiv: {
        name:      "plotdiv",
        title:     "R plot window",
        position:  {my: "right top", at: "left-10px top", of: "#command", collision: "fitflip"},
        resizable: true,
        width:     550,
        height:    570,
        minWidth:  200,
        minHeight: 200,
        closable: true
    }
}
function scanReset() {
    var toscan = getKeys(settings);
    for (var i = 0; i < toscan.length; i++) {
        var dialogobjs = getKeys(settings[toscan[i]]);
        if (any(dialogobjs, " == \"reset\"")) {
            settings[toscan[i]].reset(toscan[i]);
        }
    }
}
scanReset();
var pingobj = 0;
var pingstarted = false;
Shiny.addCustomMessageHandler("okquit",
    function(object) {
        close();
    }
);
Shiny.addCustomMessageHandler("ping",
    function(object) {
    }
);
Shiny.addCustomMessageHandler("packages",
    function(object) {
        packages = copy(object);
        refresh_cols("load");
    }
);
Shiny.addCustomMessageHandler("tempdirfile",
    function(object) {
        responseR = true;
        tempdirfile = copy(object);
    }
);
Shiny.addCustomMessageHandler("tempdatainfo",
    function(object) {
        responseR = true;
        tempdatainfo = copy(object);
    }
);
Shiny.addCustomMessageHandler("dirfile",
    function(object) {
        responseR = true;
        dirfile = copy(object);
    }
);
Shiny.addCustomMessageHandler("fullinfo",
    function(object) {
        responseR = true;
            info = object.infobjs;
            if ($("#minimize").length) {
                refresh_cols("minimize");
                checkeqtt();
            }
        if (commobj.read_table.objname !== "") {
            click_col("data_editor", "dataset", commobj.read_table.objname, others = false);
        }
        scrolleftop[commobj["data_editor"].dataset] = {"vertical": 0, "horizontal": 0};
        print_data();
        outres.console = copy(object.console);
    }
);
Shiny.addCustomMessageHandler("scrollData", 
    function(object) {
        var datasets = getKeys(object);
        for (var i = 0; i < datasets.length; i++) {
            info["data"][datasets[i]].theData = copy(object[datasets[i]].theData);
            info["data"][datasets[i]].scrollvh = copy(object[datasets[i]].scrollvh);
            info["data"][datasets[i]].dataCoords = copy(object[datasets[i]].dataCoords);
        }
        responseR = true;
    }
);
Shiny.addCustomMessageHandler("calibrate",
    function(object) {
        responseR = true;
        outres = copy(object);
    }
);
Shiny.addCustomMessageHandler("fuzzyvals",
    function(object) {
        poinths.fuzzyvals = copy(object);
        drawPointsAndThresholds();
    }
);
Shiny.addCustomMessageHandler("dataPoints",
    function(object) {
        responseR = true;
        poinths = copy(object);
    }
);
Shiny.addCustomMessageHandler("xyplot",
    function(object) {
        responseR = true;
        xyplotdata = copy(object);
    }
);
Shiny.addCustomMessageHandler("Rcommand",
    function(object) {
        if ($("#tempdiv").length) {
            input.remove();
            $("#tempdiv").remove();
        }
        var tt = extractt(object);
        if (getKeys(tt).length > 0) {
            ttfromR = tt;
            if ($("#venn").length) {
                papers["venn"]["main"].customtext = "";
                draw_venn(papers["venn"]["main"]);
            }
        }
        if ($("#xyplot").length) {
            if (object.hasOwnProperty("xyplot")) {
                xyplotdata = copy(object.xyplot);
                papers["xyplot"]["main"].xyplotdata = copy(object.xyplot);
                scaleplot(papers["xyplot"]["main"]) 
                createLabels(papers["xyplot"]["main"]);
            }
        }
        var toadd, tomodify, keys, dataset;
        var objtype = ["data", "tt", "qmc"];
        var refresh = false;
        if (object.modified !== null) {
            refresh = true;
            for (var i = 0; i < 3; i++) {
                if (object.modified[objtype[i]] !== null) {
                    if (info[objtype[i]] === null) {
                        info[objtype[i]] = {};
                    }
                    tomodify = getKeys(object.modified[objtype[i]]);
                    for (var j = 0; j < tomodify.length; j++) {
                        info[objtype[i]][tomodify[j]] = object.modified[objtype[i]][tomodify[j]];
                        if (objtype[i] == "data" && commobj["data_editor"].dataset == tomodify[j]) {
                            update_data();
                        }
                    }
                }
            }
        }
        if (object.added !== null) {
            refresh = true;
            for (var i = 0; i < 3; i++) {
                if (object.added[objtype[i]] !== null) {
                    if (info[objtype[i]] === null) {
                        info[objtype[i]] = {};
                    }
                    toadd = getKeys(object.added[objtype[i]]);
                    for (var j = 0; j < toadd.length; j++) {
                        info[objtype[i]][toadd[j]] = object.added[objtype[i]][toadd[j]];
                        if (objtype[i] == "data") {
                            click_col("data_editor", "dataset", toadd[j], others = false);
                            $("#data_editor_body").scrollTop(0);
                            $("#data_editor_body").scrollLeft(0);
                            $("#data_editor_rownames").scrollTop(0);
                            $("#data_editor_colnames").scrollLeft(0);
                            update_data();
                        }
                    }
                }
            }
        }
        var allcommobj = getKeys(commobj);
        var alldialogs = getKeys(colclicks);
        if (object.deleted !== null) {
            refresh = true;
            var position;
            for (var i = 0; i < 3; i++) {
                if (info[objtype[i]] !== null) {
                    keys = getKeys(info[objtype[i]]);
                    if (keys.length == 1) {
                        for (var d = 0; d < object.deleted.length; d++) {
                            if (keys.indexOf(object.deleted[d]) >= 0) {
                                info[objtype[i]] = null;
                                if (objtype[i] == "data") {
                                    commobj["data_editor"].dataset = "";
                                }
                            }
                        }
                    }
                    else {
                        info[objtype[i]] = copy(info[objtype[i]], exclude = object.deleted);
                        if (objtype[i] == "data" && object.deleted.indexOf(commobj["data_editor"].dataset) >= 0) {
                            commobj["data_editor"].dataset = "";
                        }
                    }
                }
                if (objtype[i] == "data" && info[objtype[i]] !== null && commobj["data_editor"].dataset == "") {
                    dataset = sortArray(getKeys(info[objtype[i]]))[0];
                    commobj["data_editor"].dataset = dataset;
                    click_col("data_editor", "dataset", dataset, others = false);
                }
                for (var d = 0; d < allcommobj.length; d++) {
                    if (commobj[allcommobj[d]].dataset !== "") {
                        for (var j = 0; j < object.deleted.length; j++) {
                            if (commobj[allcommobj[d]].dataset == object.deleted[j]) {
                                commobj[allcommobj[d]].dataset = "";
                                if (allcommobj[d] == "calibrate") {
                                    thsetter_content.remove();
                                }
                            }
                        }
                    }
                }
                for (d = 0; d < alldialogs.length; d++) {
                    if (colclicks[alldialogs[d]].dataset !== void 0) {
                        if (colclicks[alldialogs[d]].dataset.length == 1) {
                            if (object.deleted.indexOf(getKeys(colclicks[alldialogs[d]].dataset)) >= 0) {
                                colclicks = copy(colclicks, exclude = [alldialogs[d]]);
                            }
                        }
                        else {
                            if (any(colclicks[alldialogs[d]].dataset, "== true")) {
                                if (any(object.deleted, "== \"" + getTrueKeys(colclicks[alldialogs[d]].dataset)[0] + "\"")) {
                                    var allkeys = getKeys(colclicks[alldialogs[d]]);
                                    colclicks[alldialogs[d]] = copy(colclicks[alldialogs[d]], exclude = copy(allkeys, exclude = ["dataset"]));
                                }
                            }
                            colclicks[alldialogs[d]].dataset = copy(colclicks[alldialogs[d]].dataset, exclude = object.deleted);
                        }
                    }
                }
            }
        }
        if (refresh) {
            refresh_cols("all");
        }
        if (object.hasOwnProperty("poinths") && $("#calibrate").length) {
            if (commobj["calibrate"].dataset == object.poinths.dataset && commobj["calibrate"].x == object.poinths.x) {
                poinths.vals = copy(object.poinths.vals);
                poinths.prettyx = copy(object.poinths.prettyx);
                if (object.poinths.thvals !== void 0) {
                    poinths.thvals = copy(object.poinths.thvals);
                }
                drawPointsAndThresholds();
            }
        }
        if ($("#minimize").length) {
            checkeqtt();
        }
        if ($("#findRows").length) {
            refresh_cols("findRows");
        }
        if ($("#data_editor").length) {
            print_data(); 
        }
        if (object.plot) {
            plotopen = true;
            if ($("#plotdiv").length) {
                $("#plotdiv").show();
                $("#plotdiv").appendTo(document.body);
                $("#result").appendTo(document.body);
            }
            else {
                createDialog(settings["plotdiv"]);
            }
            resizePlot();
            if ($("#saveRplot").length) {
                $("#saveRplot_preview").css({
                    "background-image": "url('css/images/plot.svg?" + new Date().getTime() + "')", 
                    "background-size": "100% 100%"
                });
            }
        }
        if (history[histindex - 1] != string_command) {
            history[histindex] = string_command;
        }
        histindex = history.length;
        outres = copy(object, exclude = ["xyplot"]);
        responseR = true;
    }
);
Shiny.addCustomMessageHandler("getChanges",
    function(object) {
        responseR = true;
        outres = copy(object);
    }
);
Shiny.addCustomMessageHandler("resizePlot",
    function(object) {
        responseR = true;
        outres = copy(object);
    }
);
var shiftkey;
$(document).keydown(function (evt) {
    evt = evt || event; 
    if (evt.keyCode == 16) {
        shiftkey = true;
    }
});
$(document).keyup(function (evt) {
    evt = evt || event; 
    if (evt.keyCode == 16) {
        shiftkey = false;
    }
});
$("body").on("focus", "input, textarea", function() {
    $(this).on('keyup', function(evt) {
        evt = evt || event; 
        if (evt.keyCode == 27) { 
            tastaRcommand = "escape";
            tasta = "escape";
            input.blur();
        }
    });
    $(this).on('keydown', function(evt) {
        evt = evt || event; 
        var key = evt.which || evt.keyCode;
        if (key == 13) { 
            tastaRcommand = "enter";
            tasta = "enter";
            evt.preventDefault();
            input.blur();
        }
    });
});
function extractt(obj) {
    var toreturn = {}, keys;
    if (obj.added !== null) {
        if (obj.added.tt !== null) {
            keys = getKeys(obj.added.tt);
            toreturn = obj.added.tt[keys[0]];
        }
        else if (obj.added.qmc !== null) {
            keys = getKeys(obj.added.qmc);
            toreturn = obj.added.qmc[keys[0]];
        }
    }
    else if (obj.modified !== null) {
        if (obj.modified.tt !== null) {
            keys = getKeys(obj.modified.tt);
            toreturn = obj.modified.tt[keys[0]];
        }
        else if (obj.modified.qmc !== null) {
            keys = getKeys(obj.modified.qmc);
            toreturn = obj.modified.qmc[keys[0]];
        }
    }
    return(toreturn);
}
function talkToR() {
    reset_outres();
    Rcommand.counter = 1 - Rcommand.counter;
    Rcommand.thinfo = thinfo; 
    if (info["data"] !== null) {
        var datasets = getKeys(info["data"]);
        for (var i = 0; i < datasets.length; i++) {
            Rcommand.scrollvh[datasets[i]] = info["data"][datasets[i]].scrollvh;
        }
    }
    Rcommand.calibrate = {};
    Rcommand.xyplot = {};
    if ($("#calibrate").length) {
        if (commobj.calibrate.dataset != "" && commobj.calibrate.x != "") {
            Rcommand.calibrate = commobj.calibrate;
        }
    }
    if ($("#xyplot").length) {
        if (commobj.xyplot.dataset != "" && commobj.xyplot.x != "" && commobj.xyplot.y != "") {
            Rcommand.xyplot = commobj.xyplot;
        }
    }
    responseR = false;
    Shiny.onInputChange("Rcommand", Rcommand);
    updatecounter = 0;
    printRcommand();
}
function createCommandPromptInRconsole(prompt) {
    if (prompt === undefined) {
        prompt = ">"
    }
    var color = "#932192";
    if (prompt == "+") {
        color = "blue";
    }
    var tempdiv = "<div id = 'tempdiv'><table width = '" + ($("#result_main").width() - 8) + "'><tr><td valign = 'top' width = '15px'><span style='color:" + color + "; font-size: 14px'>" + prompt + " </span></td><td valign = 'top'><div id='txtcommand' style = 'height:19px; color:blue; font-size: 14px; padding: 2px'>" + txtcommand + "</div></td></tr></table></div>"
    $("#result_main").append(tempdiv);
    var position = $("#txtarea").position();
    $("#tempdiv").click(function(event) {
        $("#txtcommand").html("");
        showDialogToFront(settings["result"]);
        event.stopPropagation();
        input = document.createElement("textarea");
        input.value = txtcommand;
        input.setAttribute("style", getStyle());
        input.id = "txtarea";
        $("#txtcommand").append(input);
        $("#txtarea").height(getTextHeight(txtcommand, cpl.result) - 1);
        $("#result_main").scrollTop($("#result_main")[0].scrollHeight);
        var crpos;
        $("#txtarea").on("keydown", function(evt) {
            crpos = caretPosition($("#txtarea"));
        });
        $("#txtarea").on("keyup", function(evt) {
            if (evt.keyCode == 38) { 
                if (histindex > 0 && (crpos == 0 || crpos == input.value.length)) { 
                    histindex -= 1;
                    input.value = history[histindex];
                }
            }
            else if (evt.keyCode == 40) { 
                if (crpos == input.value.length) {
                    if (histindex < history.length - 1) {
                        histindex += 1;
                        input.value = history[histindex];
                    }
                    else {
                        histindex = history.length;
                        input.value = activecommand;
                    }
                }
            }
            else if (evt.keyCode < 37 || evt.keyCode > 40) {
                activecommand = input.value;
            }
            var tempheight = getTextHeight(input.value, cpl.result) - 1;
            $("#tempdiv").height(tempheight);
            $("#txtarea").height(tempheight);
            $("#result_main").scrollTop($("#result_main")[0].scrollHeight);
        });
        $("#txtarea").click(function(event) {
            event.stopPropagation();
        });
        tastaRcommand = "";
        input.addEventListener("blur", function(e) {
            if (tastaRcommand == "enter") {
                pingit();
                txtcommand = "";
                if (input.value.slice(0, 1) == "\n") {
                    input.value = input.value.slice(1);
                }
                var inputext = input.value.replace(/\r\n?/g, '\n');
                if (inputext == "" && tempcommand == "") {
                    $("#tempdiv").remove();
                    input.remove();
                    $("#result_main").append("<span style='color:#932192'>" + ((tempcommand == "")?">":"+") + " </span><br><br>");
                    tempcommand = "";
                    createCommandPromptInRconsole(prompt); 
                    $("#tempdiv").click();
                }
                else {
                    txtcommand = "";
                    string_command = inputext;
                    Rcommand.command = tempcommand + inputext;
                    talkToR();
                }
            }
            else if (tastaRcommand == "escape") {
                input.remove();
                $("#tempdiv").remove();
                if (tempcommand != "") {
                    prompt = "+";
                    color = "blue";
                }
                $("#result_main").append("<span style='color:" + color + "'>" + prompt + " </span><br><br>");
                tempcommand = "";
                createCommandPromptInRconsole();
                $("#tempdiv").click();
            }
            else {
                txtcommand = input.value;
                input.remove();
                $("#tempdiv").remove();
                createCommandPromptInRconsole(prompt);
                $("#result_main").scrollTop($("#result_main")[0].scrollHeight);
            }
            activecommand = "";
            tastaRcommand = "";
        });
        input.focus();
    });
}
function adjustDataEditorSize() {
    var vscrollbar = false, hscrollbar = false;
    var dataset = commobj["data_editor"].dataset;
    if (dataset != "" && info["data"] !== null) {
        vscrollbar = info["data"][dataset].nrows - 1 > visiblerows;
        hscrollbar = info["data"][dataset].ncols - 1 > visiblecols;
    }
    $("#data_editor").width( (visiblecols + 2)*70 + 1*(vscrollbar?scrollbarsWH:0));
    $("#data_editor").height((visiblerows + 3)*20 + 1*(hscrollbar?scrollbarsWH:0));
    $("#data_editor_main").width($("#data_editor").width());
    $("#data_editor_main").height($("#data_editor").height() - 20);
    $("#data_editor_colnames").width((visiblecols + 1)*70);
    $("#data_editor_rownames").height((visiblerows + 1)*20);
    $("#data_editor_body").width((visiblecols + 1)*70 + 1*(vscrollbar?scrollbarsWH:0));
    $("#data_editor_body").height((visiblerows + 1)*20 + 1*(hscrollbar?scrollbarsWH:0));
    if (dataset != "" && info["data"] !== null) {
        $(papers["data_editor"]["body"].canvas).width(70*info["data"][dataset].ncols);
        $(papers["data_editor"]["body"].canvas).height(20*info["data"][dataset].nrows);
        $(papers["data_editor"]["rownames"].canvas).height(20*info["data"][dataset].nrows);
        $(papers["data_editor"]["colnames"].canvas).width(70*info["data"][dataset].ncols);
    }
}
function console_command(type) {
    current_command = type;
    string_command = "";
    objname = "";
    if (type == "import") {
        if (dirfile.filepath != "") {
            if (commobj.read_table.objname != "") {
                string_command = commobj.read_table.objname + " <- ";
            }
            if (commobj.read_table.sep == ",") {
                string_command = string_command + "read.csv(\"" + 
                dirfile.filepath[0][0].replace(/\s/g, "≠") + "\"" + 
                (commobj.read_table.header?"":", header = FALSE") + 
                ((commobj.read_table.dec == ",")?", dec = \",\"":"");
            }
            else {
                string_command = string_command + "read.table(\"" + 
                dirfile.filepath[0][0].replace(/\s/g, "≠") + "\", sep = \"" +
                ((commobj.read_table.sep == "tab")?"\\t":commobj.read_table.sep) + "\"" +
                (commobj.read_table.header?", header = TRUE":"") + 
                ((commobj.read_table.dec == ",")?", dec = \",\"":"");
            }
            if (commobj.read_table.row_names.length) {
                string_command = string_command +
                ", row.names = " + (
                    (commobj.read_table.row_names % 2 >= 0)?(commobj.read_table.row_names):("\"" + commobj.read_table.row_names + "\"")
                )
            }
            string_command = string_command + ", fileEncoding = \"" + commobj.read_table.fileEncoding + "\")"
        }
    }
    else if (type == "load") {
        if (commobj["load"]["pkgdata"] != "") {
            string_command = "data(" + commobj["load"]["pkgdata"] + ")";
        }
    }
    else {
        if (info["data"] !== null) {
            if (type == "export" && commobj["export"].filename != "") {
                string_command = "export(" + commobj.export.dataset + ", file = ";
                if (commobj["export"].newfile) {
                    string_command += "\"" + (dirfile.wd + "/" + commobj["export"].filename).replace(/\s/g, "≠") + "\"";
                }
                else {
                    string_command += "\"" + (dirfile.wd + "/" + commobj["export"].filename).replace(/\s/g, "≠") + "\"";
                }
                if (commobj["export"].sep != ",") {
                    string_command += ", sep = \"" + ((commobj["export"].sep == "tab")?"\\t":commobj["export"].sep) + "\"";
                }
                if (!commobj["export"].header) {
                    string_command += ", col.names = FALSE";
                }
                else {
                    if (commobj["export"].caseid != "cases") {
                        string_command += ", caseid = \"" + commobj["export"].caseid + "\"";
                    }
                }
                string_command += ")"
            }
            if (type == "tt") {
                var outcome = getTrueKeys(colclicks.tt.outcome);
                var conditions = getTrueKeys(colclicks.tt.conditions);
                if (commobj.tt.nameit && commobj.tt.objname !== "") {
                    string_command += commobj.tt.objname + " <- ";
                }
                string_command += "truthTable(" + commobj.tt.dataset;
                if (outcome.length > 0) {
                    string_command += ", outcome = \"";
                    if (commobj.tt.neg_out) {
                        string_command += "~";
                    }
                    string_command += outcome + "\"";
                }
                if (conditions.length > 0) {
                    string_command += ", conditions = \""
                    for (var i = 0; i < conditions.length; i++) {
                        string_command += conditions[i] + ((i == conditions.length - 1)?"\"":", ");
                    }
                }
                if (commobj.tt.ic1 != "1" || (commobj.tt.ic0 != "1" && commobj.tt.ic0 != "")) {
                    if (commobj.tt.ic0 == "") {
                        string_command += ", incl.cut = " + commobj.tt.ic1;
                    }
                    else {
                        string_command += ", incl.cut = \"" + commobj.tt.ic1 + ", " + commobj.tt.ic0 + "\"";
                    }
                }
                if (commobj.tt.n_cut != "1") {
                    string_command += ", n.cut = " + commobj.tt.n_cut;
                }
                if (commobj.tt.pri_cut != "") {
                    string_command += ", pri.cut = " + commobj.tt.pri_cut;
                }
                if (commobj.tt.exclude != "") {
                    string_command += ", exclude = " + commobj.tt.exclude;
                }
                if (commobj.tt.complete) {
                    string_command += ", complete = TRUE";
                }
                if (commobj.tt.use_letters) {
                    string_command += ", use.letters = TRUE";
                }
                if (commobj.tt.show_cases) {
                    string_command += ", show.cases = TRUE";
                }
                if (commobj.tt.show_cases && commobj.tt.dcc) {
                    string_command += ", dcc = TRUE";
                }
                var sorts = getTrueKeys(commobj.tt.sort_sel);
                if (sorts.length > 0) {
                    string_command += ", sort.by = \"";
                    for (var i = 0; i < sorts.length; i++) {
                        string_command += sorts[i] +
                                          ((commobj.tt.sort_by[sorts[i]])?"":"+") +
                                          ((i == sorts.length - 1)?"\"":", ");
                    }
                }
                string_command += ")";
            }
            if (type == "calibrate") {
                var col = (getKeys(colclicks).indexOf("calibrate") >= 0)?getTrueKeys(colclicks.calibrate.x):"";
                if (col.length > 0) { 
                    if (!commobj.calibrate.same && commobj.calibrate.newvar != "") {
                        string_command = commobj.calibrate.dataset + "$" + commobj.calibrate.newvar;
                    }
                    else {
                        string_command = commobj.calibrate.dataset + "$" + col;
                    }
                    string_command += " <- calibrate(" + commobj.calibrate.dataset + "$" + col;
                    if (commobj.calibrate.type == "crisp") {
                        string_command += ", type = \"crisp\"";
                    }
                    if (commobj.calibrate.thresholds.length > 0) {
                        var valid = new Array();
                        for (var i = 0; i < commobj.calibrate.thresholds.length; i++) {
                            if (commobj.calibrate.thresholds[i] != "" && commobj.calibrate.thresholds[i] != void 0) {
                                valid.push(i);
                            }
                        }
                        if (valid.length > 0) {
                            if (valid.length == 1) {
                                if (commobj.calibrate.type == "crisp") {
                                   string_command += ", thresholds = " + commobj.calibrate.thresholds[valid[0]];
                                }
                            }
                            else {
                                if (commobj.calibrate.type == "crisp") {
                                    string_command += ", thresholds = c(";
                                    for (var i = 0; i < valid.length; i++) {
                                        string_command += commobj.calibrate.thresholds[valid[i]] + ((i < valid.length - 1)?", ":")");
                                    }
                                }
                                else {
                                    if ((valid.length == 3 && commobj.calibrate.thnames[0].substring(1, 2) != "1") || valid.length == 6) {
                                        string_command += ", thresholds = \"";
                                        for (var i = 0; i < valid.length; i++) {
                                            string_command += (commobj.calibrate.thnames[valid[i]] + "=") + commobj.calibrate.thresholds[valid[i]] + ((i == valid.length - 1)?"\"":", ");
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (commobj.calibrate.type == "crisp") {
                        if (commobj.calibrate.include) {
                            string_command += ", include = TRUE";
                        }
                    }
                    if (commobj.calibrate.type == "fuzzy") {
                        if (commobj.calibrate.end) {
                            if (!commobj.calibrate.logistic) {
                                string_command += ", logistic = FALSE";
                                if (commobj.calibrate.ecdf) {
                                    string_command += ", ecdf = TRUE";
                                }
                            }
                            else if (commobj.calibrate.idm != "0.95") {
                                string_command += ", idm = " + commobj.calibrate.idm;
                            }
                        }
                    }
                    if (commobj.calibrate.below != "1") {
                        string_command += ", below = " + commobj.calibrate.below;
                    }
                    if (commobj.calibrate.above != "1") {
                        string_command += ", above = " + commobj.calibrate.above;
                    }
                    string_command += ")";
                }
                else {
                    string_command = "";
                }
            }
            if (type == "recode") {
                var col = (getKeys(colclicks).indexOf("recode") >= 0)?getTrueKeys(colclicks.recode.x):"";
                var uniques = unique(commobj.recode.newv);
                if (col.length > 0) {
                    if (!commobj.recode.same && commobj.recode.newvar != "") {
                        string_command = commobj.recode.dataset + "$" + commobj.recode.newvar;
                    }
                    else {
                        string_command = commobj.recode.dataset + "$" + col;
                    }
                    string_command += " <- recode(" + commobj.recode.dataset + "$" + col + ", \"";
                    var nl = commobj.recode.newv.length; 
                    var temp, oldvals;
                    for (var i = 0; i < uniques.length; i++) {
                        temp = new Array();
                        oldvals = "";
                        for (var j = 0; j < nl; j++) {
                            if (commobj.recode.newv[j] == uniques[i]) {
                                temp.push(commobj.recode.oldv[j]);
                            }
                        }
                        for (j = 0; j < temp.length; j++) {
                            oldvals += temp[j] + ((j < temp.length - 1)?",":"")
                        }
                        string_command += oldvals + "=" + uniques[i] + ((i < uniques.length - 1)?"; ":"");
                    }
                    string_command += "\")";
                }
                else {
                    string_command = "";
                }
            }
        }
        if (type == "findRows") {
            var countchecked = 0;
            for (var i = 0; i < 4; i++) {
                countchecked += commobj.findRows.type[i];
            }
            if (commobj.findRows.obj != "" && countchecked > 0) {
                if (commobj.findRows.objname !== "") {
                    string_command += commobj.findRows.objname + " <- ";
                }
                string_command += "findRows(";
                if (commobj.findRows.type[0] || commobj.findRows.type[1]) {
                    if (commobj.findRows.expression != "") {
                        string_command += ", expression = \"" + commobj.findRows.expression + "\", " + commobj.findRows.obj;
                    }
                }
                else {
                    string_command += "obj = " + commobj.findRows.obj;
                }
                if (!commobj.findRows.remainders) {
                    string_command += ", remainders = FALSE";
                }
                if (commobj.findRows.type[0]) {
                    string_command += ", type = 0";
                }
                else {
                    countchecked = 0;
                    for (var i = 1; i < 4; i++) {
                        countchecked += commobj.findRows.type[i];
                    }
                    if (countchecked == 1) {
                        for (var i = 2; i < 4; i++) {
                            if (commobj.findRows.type[i]) {
                                string_command += ", type = ";
                                string_command += i;
                            }
                        }
                    }
                    else {
                        string_command += ", type = c(";
                        var first = true;
                        for (var i = 1; i < 4; i++) {
                            if (commobj.findRows.type[i]) {
                                if (first) {
                                    string_command += i;
                                    first = false;
                                }
                                else {
                                    string_command += ", " + i;
                                }
                            }
                        }
                        string_command += ")";
                    }
                }
                string_command += ")";
            }
        }
        if (info["data"] !== null || info["tt"] !== null) {
            if (type == "minimize") {
                var outcome = new Array();
                var conditions = new Array();
                if (colclicks.minimize !== void 0) {
                    if (colclicks.minimize.outcome !== void 0) {
                        outcome = getTrueKeys(colclicks.minimize.outcome);
                    }
                    if (colclicks.minimize.conditions !== void 0) {
                        conditions = getTrueKeys(colclicks.minimize.conditions);
                    }
                }
                if (commobj.minimize.objname !== "") {
                    string_command += commobj.minimize.objname + " <- ";
                }
                string_command += "minimize(" + commobj.minimize.dataset;
                if (outcome.length > 0 && commobj.minimize.source == "data") {
                    string_command += ", outcome = \"";
                    if (commobj.minimize.neg_out && outcome.length == 1) {
                        string_command += "~";
                    }
                    for (var i = 0; i < outcome.length; i++) {
                        string_command += outcome[i] + ((i == outcome.length - 1)?"\"":", ");
                    }
                }
                if (conditions.length > 0 && commobj.minimize.source == "data") {
                    string_command += ", conditions = \"";
                    for (var i = 0; i < conditions.length; i++) {
                        string_command += conditions[i] + ((i == conditions.length - 1)?"\"":", ");
                    }
                }
                if (commobj.minimize.use_letters) {
                    string_command += ", use.letters = TRUE";
                }
                if (commobj.minimize.ic1 != "1") {
                    if (commobj.minimize.ic0 == "") {
                        string_command += ", incl.cut = " + commobj.minimize.ic1;
                    }
                    else {
                        string_command += ", incl.cut = \"" + commobj.minimize.ic1 + ", " + commobj.minimize.ic0 + "\"";
                    }
                }
                if (commobj.minimize.pri_cut != "") {
                    string_command += ", pri.cut = " + commobj.minimize.pri_cut;
                }
                if (commobj.minimize.n_cut != "1") {
                    string_command += ", n.cut = " + commobj.minimize.n_cut;
                }
                if (commobj.minimize.include[0] != "" || commobj.minimize.include[1] != "") {
                    string_command += ", include = \"";
                    string_command += commobj.minimize.include[0] + (
                        (commobj.minimize.include[1] != "") ? (
                            ((commobj.minimize.include[0] != "") ? ", " : "") +
                        commobj.minimize.include[1]) : "");
                    string_command += "\"";
                }
                if (commobj.minimize.dir_exp != "" && commobj.minimize.include[0] == "?") {
                    string_command += ", dir.exp = \"" + commobj.minimize.dir_exp + "\"";
                }
                if (commobj.minimize.details) {
                    string_command += ", details = TRUE";
                }
                if (commobj.minimize.all_sol) {
                    string_command += ", all.sol = TRUE";
                }
                if (commobj.minimize.first_min) {
                    string_command += ", first.min = TRUE";
                }
                if (commobj.minimize.row_dom) {
                    string_command += ", row.dom = TRUE";
                }
                if (commobj.minimize.pi_cons != "0") {
                    string_command += ", pi.cons = " + commobj.minimize.pi_cons;
                }
                if (commobj.minimize.pi_depth != "0") {
                    string_command += ", pi.depth = " + commobj.minimize.pi_depth;
                }
                if (commobj.minimize.sol_cons != "0") {
                    string_command += ", sol.cons = " + commobj.minimize.sol_cons;
                }
                if (commobj.minimize.sol_cov != "1") {
                    string_command += ", sol.cov = " + commobj.minimize.sol_cov;
                }
                if (commobj.minimize.sol_depth != "0") {
                    string_command += ", sol.depth = " + commobj.minimize.sol_depth;
                }
                string_command += ")";
            }
        }
    }
    $("#command_main").html(strwrap(string_command, cpl.command));
}
function click_col(dialog, identifier, col, others) {
    if (missing(others)) {
        others = true;
    }
    if (getKeys(colclicks).indexOf(dialog) < 0) {
        colclicks[dialog] = new Object;
    }
    if (getKeys(colclicks[dialog]).indexOf(identifier) < 0) {
        colclicks[dialog][identifier] = new Object;
    }
    else {
        if (!others) {
            others = getKeys(colclicks[dialog][identifier]);
            for (var i = 0; i < others.length; i++) {
                colclicks[dialog][identifier][others[i]] = false;
            }
        }
    }
    colclicks[dialog][identifier][col] = true;
    commobj[dialog][identifier] = col;
}
function checkeqtt() {
    if ($("#minimize").length) {
        if (commobj.minimize.source == "data") {
            papers["minimize"]["main"].ic1.attr({"text": commobj.minimize.ic1});
            papers["minimize"]["main"].ic0.attr({"text": commobj.minimize.ic0});
            papers["minimize"]["main"].pri_cut.attr({"text": commobj.minimize.pri_cut});
            papers["minimize"]["main"].n_cut.attr({"text": commobj.minimize.n_cut});
            papers["minimize"]["main"].neg_out.refresh(commobj.minimize.neg_out);
            papers["minimize"]["main"].show_cases.refresh(commobj.minimize.show_cases);
            papers["minimize"]["main"].use_letters.refresh(commobj.minimize.use_letters);
        }
        else if (commobj.minimize.source == "tt" && commobj.minimize.dataset != "") {
            if (info["tt"][commobj.minimize.dataset] !== void 0) {
                papers["minimize"]["main"].ic1.attr({"text": info["tt"][commobj.minimize.dataset].options["incl.cut"][0]});
                papers["minimize"]["main"].ic0.attr({"text": ""});
                if (info["tt"][commobj.minimize.dataset].options["incl.cut"].length == 2) {
                    papers["minimize"]["main"].ic0.attr({"text": info["tt"][commobj.minimize.dataset].options["incl.cut"][1]});
                }
                papers["minimize"]["main"].pri_cut.attr({"text": info["tt"][commobj.minimize.dataset].options["pri.cut"]});
                papers["minimize"]["main"].n_cut.attr({"text": info["tt"][commobj.minimize.dataset].options["n.cut"]});
                papers["minimize"]["main"].neg_out.refresh(info["tt"][commobj.minimize.dataset].options["neg.out"]);
                papers["minimize"]["main"].show_cases.refresh(info["tt"][commobj.minimize.dataset].options["show.cases"]);
                papers["minimize"]["main"].use_letters.refresh(info["tt"][commobj.minimize.dataset].options["use.letters"]);
            }
        }
    }
}
function print_cols(dialog, identifier, options) {
    var selection = options.selection;
    var cols = options.cols;
    var selectable = options.selectable;
    var objtype = (dialog == "minimize")?(commobj["minimize"].source):"data";
    var paper = papers[dialog][identifier];
    paper.clear();
    var numerics, calibrated;
    if (cols.length > 0) {
        if (["package", "pkgdata", "dataset", "tt", "qmc"].indexOf(identifier) >= 0) {
            cols = sortArray(options.cols);
        }
        if (dialog !== "load" && dialog !== "findRows" && commobj[dialog] !== void 0) {
            if (commobj[dialog].dataset != "") {
                numerics = info[objtype][commobj[dialog].dataset].numerics;
                calibrated = info[objtype][commobj[dialog].dataset].calibrated; 
            }
        }
        if (selection === void 0) {
            selection = "single";
        }
        if (getKeys(colclicks).indexOf(dialog) < 0) {
            colclicks[dialog] = new Object;
        }
        if (getKeys(colclicks[dialog]).indexOf(identifier) < 0) {
            colclicks[dialog][identifier] = new Object;
            for (var i = 0; i < cols.length; i++) {
                colclicks[dialog][identifier][cols[i]] = false;
            }
        }
        else {
            for (var i = 0; i < cols.length; i++) {
                if (getKeys(colclicks[dialog][identifier]).indexOf(cols[i]) < 0) {
                    colclicks[dialog][identifier][cols[i]] = false;
                }
            }
        }
        if (dialog == "minimize" && commobj["minimize"].source == "tt") {
            if (identifier != "dataset") {
                var ocheck = getTrueKeys(colclicks[dialog][identifier]);
                if (ocheck.length > 0) {
                    for (var i = 0; i < ocheck.length; i++) {
                        colclicks[dialog][identifier][ocheck[i]] = false;
                    }
                }
            }
            if (identifier == "outcome") {
                var outcome = info["tt"][commobj[dialog].dataset].options.outcome;
                colclicks[dialog][identifier][outcome] = true;
                commobj[dialog].outcome = [outcome];
            }
            else if (identifier == "conditions") {
                var conditions = info["tt"][commobj[dialog].dataset].options.conditions;
                for (var i = 0; i < conditions.length; i++) {
                    colclicks[dialog][identifier][conditions[i]] = true;
                }
                commobj[dialog].conditions = conditions;
            }
        }
        canvas_height = cols.length * 20;
        var rects_back = new Array(cols.length);
        var texts = new Array(cols.length);
        var rects = new Array(cols.length);
        var clicks = [-1, -1];
        if (selection == "none") {
            var colset = paper.set();
            if (Array.isArray(cols)) {
                for (var i = 0; i < cols.length; i++) {
                    colset.push(paper.rect(0, i*20 + 0.5, 220, 19).attr({fill: colclicks[dialog][identifier][cols[i]]?"#79a74c":"#ffffff", stroke: "none"}));
                    colset.push(paper.text(10, 10 + i*20, cols[i]).attr({"text-anchor": "start", "font-size": "14px", fill: colclicks[dialog][identifier][cols[i]]?"white":"black"}));
                }
            }
            else {
                colset.push(sat(paper.text(10, 11, cols)));
            }
            var colsetbox = colset.getBBox();
            $(paper.canvas).width(colsetbox.width + 20); 
            if (dialog == "minimize") {
                paper.rect(0, 0, 138, 120).attr({fill: "#000000", stroke: "none", "fill-opacity": "0.2"});
            }
        }
        else {
            if (dialog == "load" || info[objtype] !== null) {
                for (var i = 0; i < cols.length; i++) {
                    rects_back[i] = paper.rect(0, i*20 + 0.5, 220, 19).attr({fill: colclicks[dialog][identifier][cols[i]]?"#79a74c":"#ffffff", stroke: "none"});
                    texts[i] = paper.text(10, 10 + i*20, cols[i]).attr({"text-anchor": "start", "font-size": "14px", fill: colclicks[dialog][identifier][cols[i]]?"white":"black"});
                    opacity = 0;
                    if (selectable.indexOf("numerics") >= 0) {
                        if (!numerics[i]) {
                            opacity = 0.3;
                        }
                    }
                    if (selectable.indexOf("calibrated") >= 0) {
                        if (!calibrated[i]) {
                            opacity = 0.3;
                        }
                    }
                    rects[i] = paper.rect(0, i*20, 220, 20).attr({fill: "#333333", stroke: "none", "fill-opacity": opacity});
                    rects[i].selectable = (opacity == 0);
                    rects[i].id = i;
                    rects[i].name = cols[i];
                    if (rects[i].selectable) {
                        rects[i].click(function(event) {
                            if (selection == "multiple") {
                                if (event.shiftKey) {
                                    if (clicks[0] > -1) {
                                        var x1 = clicks[0];
                                        var x2 = this.id;
                                        if (x1 > x2) {
                                            x1 = this.id;
                                            x2 = clicks[0];
                                        }
                                        var firstvar = rects[clicks[0]].name;
                                        for (var k = x1; k < x2 + 1; k++) {
                                            if (rects[k].selectable) {
                                                rects_back[k].attr({fill: colclicks[dialog][identifier][firstvar]?"#79a74c":"#ffffff", stroke: "none"});
                                                texts[k].attr({"text-anchor": "start", "font-size": "14px", fill: colclicks[dialog][identifier][firstvar]?"white":"black"});
                                                colclicks[dialog][identifier][rects[k].name] = colclicks[dialog][identifier][firstvar];
                                            }
                                        }
                                    }
                                    else {
                                        clicks[0] = this.id;
                                        rects_back[this.id].attr({fill: "#79a74c", stroke: "none"});
                                        texts[this.id].attr({"text-anchor": "start", "font-size": "14px", fill: "white"});
                                        colclicks[dialog][identifier][this.name] = true;
                                        if (dialog == "recode") {
                                            checkRecodeSelections(colclicks, papers["recode"]["main"]);
                                        }
                                    }
                                }
                                else { 
                                    clicks[0] = this.id;
                                    colclicks[dialog][identifier][this.name] = !colclicks[dialog][identifier][this.name];
                                    rects_back[this.id].attr({fill: colclicks[dialog][identifier][this.name]?"#79a74c":"#ffffff", stroke: "none"});
                                    texts[this.id].attr({"text-anchor": "start", "font-size": "14px", fill: colclicks[dialog][identifier][this.name]?"white":"black"});
                                    if (dialog == "recode") {
                                        checkRecodeSelections(colclicks, papers["recode"]["main"]);
                                    }
                                }
                                console_command(dialog);
                            }
                            else if (selection == "single") {
                                if (this.selectable) {
                                    if (!colclicks[dialog][identifier][this.name]) {
                                        for (var k = 0; k < cols.length; k++) {
                                            rects_back[k].attr({fill: "#ffffff", stroke: "none"});
                                            texts[k].attr({"text-anchor": "start", "font-size": "14px", fill: "black"});
                                            colclicks[dialog][identifier][rects[k].name] = false;
                                        }
                                        colclicks[dialog][identifier][this.name] = true;
                                        rects_back[this.id].attr({fill: "#79a74c", stroke: "none"});
                                        texts[this.id].attr({"text-anchor": "start", "font-size": "14px", fill: "white"});
                                        if (dialog == "load") {
                                            var former = commobj[dialog][identifier];
                                            commobj[dialog][identifier] = this.name;
                                            papers["load"]["title"].clear();
                                            if (identifier == "package") {
                                                if (former != this.name) {
                                                    commobj.load.pkgdata = "";
                                                    colclicks["load"]["pkgdata"] = new Array();
                                                    $("#load_pkgdata").scrollTop(0);
                                                }
                                                refresh_cols(dialog);
                                            }
                                            else { 
                                                var text = packages[commobj.load.package][commobj.load.pkgdata].split(" ");
                                                var rows = 0;
                                                var rowtext = "";
                                                for (var i = 0; i < text.length; i++) {
                                                    if ((rowtext + text[i]).length > 55) {
                                                        sat(papers["load"]["title"].text(10, 10 + 20*rows, rowtext));
                                                        rowtext = text[i] + " ";
                                                        rows += 1;
                                                    }
                                                    else {
                                                        rowtext += text[i] + " ";
                                                    }
                                                }
                                                sat(papers["load"]["title"].text(10, 10 + 20*rows, rowtext));
                                            }
                                        }
                                        else if (dialog == "calibrate") {
                                            for (var i = 0; i < 6; i++) {
                                                ths[i].attr({"text": ""});
                                                commobj.calibrate.thresholds[i] = "";
                                                commobj.calibrate.thnames[i] = "";
                                            }
                                            if (identifier == "x") {
                                                commobj["calibrate"].x = this.name;
                                                thinfoToR();
                                            }
                                            else if (identifier == "dataset") {
                                                resetPoinths();
                                                drawPointsAndThresholds();
                                            }
                                        }
                                        else if (dialog == "export") {
                                            commobj.export.dataset = this.name;
                                        }
                                        else if (dialog == "findRows") {
                                            commobj.findRows.obj = this.name;
                                        }
                                        else if (dialog == "xyplot") {
                                            if (identifier !== "dataset") {
                                                commobj.xyplot[identifier] = this.name;
                                                if (commobj.xyplot.x != "" && commobj.xyplot.y != "") {
                                                    updatecounter = 0;
                                                    commobj.xyplot.counter += 1;
                                                    lastvals = xyplotdata;
                                                    responseR = false;
                                                    Shiny.onInputChange("xyplot", commobj.xyplot);
                                                    doWhenXYplotPointsAreReturned();
                                                }
                                            }
                                        }
                                        else if (dialog == "data_editor") {
                                            if (commobj["data_editor"]["dataset"] != this.name) {
                                                commobj["data_editor"]["dataset"] = this.name;
                                                print_data();
                                            }
                                        }
                                    }
                                    if (dialog == "data_editor") {
                                        $("#data_editor_dataset_border").remove();
                                    }
                                    if (identifier == "dataset" && ["minimize", "tt", "calibrate", "recode", "xyplot"].indexOf(dialog) >= 0) {
                                        if (commobj[dialog][identifier] != this.name) {
                                            var identifiers = getKeys(colclicks[dialog]);
                                            for (var i = 0; i < identifiers.length; i++) {
                                                if (identifiers[i] !== "dataset") {
                                                    colclicks[dialog][identifiers[i]] = new Object;
                                                    $("#" + dialog + "_" + identifiers[i]).scrollTop(0);
                                                }
                                            }
                                        }
                                        if (dialog == "xyplot") {
                                            xyplotdata = new Array();
                                            commobj.xyplot.x = "";
                                            commobj.xyplot.y = "";
                                            draw_xyplot(papers["xyplot"]["main"]);
                                        }
                                        commobj[dialog][identifier] = this.name;
                                        refresh_cols(dialog);
                                        checkeqtt();
                                    }
                                    console_command(dialog);
                                }
                            }
                        });
                    }
                }
            }
            else {
                canvas_height = 0;
            }
        }
        $(paper.canvas).height(canvas_height);
    }
}
function strwrap(str, width, prefix, html) {
    html = void 0==html?true:html;
    prefix = void 0==prefix?"":prefix;
    str = str.replace("csv(", "£");
    str = str.replace("table(", "§");
    str = str.replace(/\s/g, "∞");
    str = str.replace(/\//g, "±");
    var toreturn;
    var crev = {"£": "csv(", "§": "table(", "∞":(html?"&nbsp;":" "), "±": "/", "≠": " "};
    if (str.length > width) {
        var q = 0;
        for (var p = width; p > 1; p--) {
            if (["£", "§", "∞", "±"].indexOf(str[p]) >= 0 && q == 0) {
                q = p;
            }
        }
        if (q > 0) {
            var left = str.substring(0, q + (str[q] != "∞"));
            var right = prefix + str.substring(q + 1);
            toreturn = left + (html?"<br>":"\n") + strwrap(right, width, prefix);
            toreturn = toreturn.replace(/£|§|∞|±|≠/g, function(x) {return crev[x]});
            return (toreturn);
        }
        else {
            var q = 0;
            for (var p = 1; p < str.length; p++) {
                if (["£", "§", "∞", "±"].indexOf(str[p]) >= 0 && q == 0) {
                    q = p;
                }
            }
            if (q > 0) {
                var left = str.substring(0, q + (str[q] != "∞"));
                var right = prefix + str.substring(q + 1);
                toreturn = left + (html?"<br>":"\n") + strwrap(right, width, prefix);
                return (toreturn.replace(/£|§|∞|±|≠/g, function(x) {return crev[x]}));
            }
            else {
                return (str.replace(/£|§|∞|±|≠/g, function(x) {return crev[x]}))
            }
        }
    }
    return (str.replace(/£|§|∞|±|≠/g, function(x) {return crev[x]}));
}
function getTextHeight(text, width) {
    let input2 = document.createElement("textarea");
    text = strwrap(text, width, "", false);
    input2.value = text;
    input2.setAttribute("style", getStyle(fresh = true));
    input2.id = "txtarea2";
    $("#command").append(input2);
    let iheight = $("#txtarea2")[0].scrollHeight;
    $("#txtarea2").remove();
    return(iheight);
}
function getTextWidth(string) {
    var paper = Raphael(0, 0, 0, 0);
        paper.canvas.style.visibility = "hidden";
    var BBox = sat(paper.text(-10, -10, string)).getBBox();
    paper.remove();
    return BBox.width;
}
function getTrimmedText(text, width) {
    var temp = "";
    var stop = false;
    for (var i = 0; i < text.length; i++) {
        if (getTextWidth(temp + text[i]) <= width) {
            temp += stop?"":text[i];
        }
        else {
            stop = true;
        }
    }
    return(temp + "...");
}
function print_data() {
    if ($("#data_editor").length) {
        adjustDataEditorSize();
        if (papers["data_editor"]["topleft"].constant === void 0) {
            papers["data_editor"]["topleft"].constant = papers["data_editor"]["topleft"].rect(0, 0, 70, 20).attr({fill: "#f2f2f2", stroke: "#d7d7d7", "fill-opacity": 1});
            papers["data_editor"]["topleft"].colsrect = -100;
            papers["data_editor"]["topleft"].colsrect_show = false;
            papers["data_editor"]["topleft"].rowsrect = -100;
            papers["data_editor"]["topleft"].rowsrect_show = false;
            papers["data_editor"]["topleft"].bodyrect = [-100, -100];
            papers["data_editor"]["topleft"].bodyrect_show = false;
            papers["data_editor"]["topleft"].constant.attr({"fill-opacity": 1, "stroke": "#d7d7d7"})
            papers["data_editor"]["topleft"].name = sat(papers["data_editor"]["topleft"].text(5, 10, "Choose"), {size: 12, "font-weight": "bold", "color": "#cb2626"});
            var choose = sat(papers["data_editor"]["topleft"].rect(0, 0, 70, 20), {"stroke": "#000", "sw": 2});
            var data = {border: true,  left:  0, top: 40, width: (70 + scrollbarsWH), height: 100};
            choose.click(function() {
                if ($("#data_editor_dataset_border").length) {
                    $("#data_editor_dataset_border").remove();
                }
                else {
                    addDiv("data_editor", "dataset", data);
                    papers["data_editor"]["dataset"] = Raphael("data_editor_dataset", data.width, data.height);
                    refresh_cols("data_editor");
                }
            });
        }
        var dataset = commobj["data_editor"].dataset;
        if (dataset != "") {
            papers["data_editor"]["body"].setSize(70*info["data"][dataset].ncols, 20*info["data"][dataset].nrows);
            papers["data_editor"]["rownames"].setSize(70                        , 20*info["data"][dataset].nrows);
            papers["data_editor"]["colnames"].setSize(70*info["data"][dataset].ncols, 20);
        }
        else {
            papers["data_editor"]["body"].setSize(70*(visiblecols + 1), 20*(visiblerows + 1));
            papers["data_editor"]["rownames"].setSize(70                  , 20*(visiblerows + 1));
            papers["data_editor"]["colnames"].setSize(70*(visiblecols + 1), 20                  );
        }
            update_data();
    }
}
function update_data() {
    if ($("#data_editor").length) {
        var activeSquare = function() {
            gridset.toFront();
            if (dataset != "") {
                datacover.toFront();
            }
            if (papers["data_editor"]["topleft"].rowsrect_show) {
                rowsrect = papers["data_editor"]["rownames"].rect(1, papers["data_editor"]["topleft"].rowsrect, 68, 18).attr({"stroke-width": 1.3});
            }
            if (papers["data_editor"]["topleft"].colsrect_show) {
                colsrect = papers["data_editor"]["colnames"].rect(papers["data_editor"]["topleft"].colsrect, 1, 68, 18).attr({"stroke-width": 1.3});
            }
            if (papers["data_editor"]["topleft"].bodyrect_show) {
                bodyrect = papers["data_editor"]["body"].rect(papers["data_editor"]["topleft"].bodyrect[0], papers["data_editor"]["topleft"].bodyrect[1], 68, 18).attr({"stroke-width": 1.3});
            }
        }
        var dataset = commobj["data_editor"].dataset;
        papers["data_editor"]["topleft"].name.attr({"text": ((dataset == "")?"Choose":dataset)});
        papers["data_editor"]["colnames"].clear();
        papers["data_editor"]["rownames"].clear();
        papers["data_editor"]["body"].clear();
        if (scrolleftop[commobj["data_editor"].dataset] !== void 0) {
            $("#data_editor_rownames").scrollTop(scrolleftop[commobj["data_editor"].dataset].vertical);
            $("#data_editor_colnames").scrollLeft(scrolleftop[commobj["data_editor"].dataset].horizontal);
            $("#data_editor_body").scrollTop(scrolleftop[commobj["data_editor"].dataset].vertical);
            $("#data_editor_body").scrollLeft(scrolleftop[commobj["data_editor"].dataset].horizontal);
        }
        var Xshift = Math.floor($("#data_editor_body").scrollLeft()/70);
        var Yshift = Math.floor($("#data_editor_body").scrollTop()/20);
        var temp, tocompare, textToPrint, tobe, temprect;
        var bodyrect = papers["data_editor"]["body"].rect(-100, 0, 70, 20);
        var colsrect = papers["data_editor"]["body"].rect(-100, 0, 70, 20);
        var rowsrect = papers["data_editor"]["body"].rect(-100, 0, 70, 20);
        var bodygridtext = "", colgridtext = "", rowgridtext = "";
        papers["data_editor"]["colnames"].rect(70*(Xshift - 25), 0, 70*(Xshift + 60), 20)
        .attr({fill: "#f2f2f2", stroke: "#d7d7d7"});
        papers["data_editor"]["rownames"].rect(0, 20*(Yshift - 50), 70, 20*(Yshift + 120))
        .attr({fill: "#f2f2f2", stroke: "#d7d7d7"});
        for (var i = Xshift - 25; i < Xshift + 60; i++) { 
            bodygridtext += "M" + 70*i + "," + 20*(Yshift - 50) + "L" + 70*i + "," + 20*(Yshift + 120);
            colgridtext += "M" + 70*i + ",0 L" + 70*i + ",20";
        }
        for (var i = Yshift - 50; i < Yshift + 120; i++) { 
            bodygridtext += "M" + 70*(Xshift - 25) + "," + 20*i + "L" + 70*(Xshift + 60) + "," + 20*i;
            papers["data_editor"]["rownames"].path("M" + 0 + "," + 20*i + "L 70" + "," + 20*i).attr({stroke: "#d7d7d7"});
        }
        gridset = papers["data_editor"]["body"].path(bodygridtext).attr({stroke: "#d7d7d7"});
        papers["data_editor"]["colnames"].path(colgridtext).attr({stroke: "#d7d7d7"});
        var getCoords = function(event) {
            var scrollX = $("#data_editor_body").scrollLeft()%70;
            var scrollY = $("#data_editor_body").scrollTop()%20;
            var mouseX = Math.floor((event.clientX + $(window).scrollLeft() - testX - 70 + scrollX)/70);
            var mouseY = Math.floor((event.clientY + $(window).scrollTop() - 3*(navigator.browserType == "Firefox") - testY - 40 + scrollY)/20);
            var Xshift = Math.floor($("#data_editor_body").scrollLeft()/70);
            var Yshift = Math.floor($("#data_editor_body").scrollTop()/20);
            var mX = mouseX - Math.round(scrollX/70);
            var mY = mouseY - Math.round(scrollY/20);
            return({
                "mX": mX,
                "mY": mY,
                "mouseX": mouseX,
                "mouseY": mouseY,
                "scrollX": scrollX,
                "scrollY": scrollY,
                "Xshift": Xshift,
                "Yshift": Yshift
            });
        }
        var scrollvh = [0, 0, visiblerows, visiblecols];
        if (info["data"] !== null && dataset != "") {
            scrollvh[0] = info["data"][dataset].scrollvh[0];
            scrollvh[1] = info["data"][dataset].scrollvh[1];
            scrollvh[2] = Math.min(visiblerows, info["data"][dataset].nrows - 1);
            scrollvh[3] = Math.min(visiblecols, info["data"][dataset].ncols - 1);
        }
        var column;
        if (dataset != "") {
            for (var i = 0; i < scrollvh[3] + 1; i++) {
                sat(papers["data_editor"]["colnames"].text(5 + 70*(i + scrollvh[1]), 10, info["data"][dataset].colnames[i + scrollvh[1]]),
                    {"clip": (70*(i + scrollvh[1])) + ", 0, 68, 20"});
            }
            for (var j = 0; j < scrollvh[2] + 1; j++) { 
                 sat(papers["data_editor"]["rownames"].text(5, 10 + 20*(j + scrollvh[0]), info["data"][dataset].rownames[j + scrollvh[0]]),
                    {"clip": "0, " + 20*(j + scrollvh[0]) + ", 68, 20"});
            }
            for (i = 0; i < scrollvh[3] + 1; i++) {
                column = info["data"][dataset].theData[i];
                for (j = 0; j < column.length; j++) {
                    if (column[j] != undefined) {
                        sat(papers["data_editor"]["body"].text(5 + 70*(i + scrollvh[1]), 10 + 20*(j + scrollvh[0]), ("" + column[j])), 
                            {"clip": 70*(i + scrollvh[1]) + ", " + 20*(j + scrollvh[0]) + ", 68, 20"});
                    }
                    else {
                        papers["data_editor"]["body"].text(5 + 70*(i + scrollvh[1]), 10 + 20*(j + scrollvh[0]), "");
                    }
                }
            }
            var colnamescover = papers["data_editor"]["colnames"].rect(0, 0, 70*info["data"][dataset].ncols, 20)
            .attr({fill: "#ffffff", stroke: "none", "fill-opacity": 0})
            .click(function(event) {
                var coords = getCoords(event);
                bodyrect.hide();
                rowsrect.hide();
                colsrect.remove();
                papers["data_editor"]["topleft"].bodyrect_show = false;
                papers["data_editor"]["topleft"].rowsrect_show = false;
                papers["data_editor"]["topleft"].colsrect_show = true;
                colsrect = papers["data_editor"]["colnames"].rect(70*(coords.mouseX + coords.Xshift) + 1, 1, 68, 18).attr({"stroke-width": 1.3});
                papers["data_editor"]["topleft"].colsrect = 70*(coords.mouseX + coords.Xshift) + 1;
                if ($("#data_editor_dataset_border").length) {
                    $("#data_editor_dataset_border").remove();
                }
            })
            .dblclick(function(event) {
                var coords = getCoords(event);
                colsrect.hide();
                temp = info["data"][dataset].colnames[coords.mouseX + coords.Xshift];
                tobe = sat(papers["data_editor"]["colnames"].text(0, 0, temp));
                papers["data_editor"]["colnames"].inlineTextEditing(tobe);
                input = tobe.inlineTextEditing.startEditing(
                    70*coords.mouseX - coords.scrollX + 70, 
                    20 - 1*(navigator.browserType == "Firefox"), 
                    70, 
                    20,
                    "from_data_editor",
                    "#f2f2f2");
                input.cover = this;
                input.addEventListener("blur", function(e) {
                    tobe.inlineTextEditing.stopEditing(tasta);
                    tocompare = tobe.attr("text");
                    if (temp != tocompare) {
                        colclicks = changeCol(colclicks, temp, tocompare);
                        refresh_cols("all", exclude = ["import"]);
                            papers["data_editor"]["colnames"].rect(70*(coords.mouseX + coords.Xshift), 0, 70, 20)
                            .attr({fill: "#f2f2f2", stroke: "#d7d7d7"});
                            sat(papers["data_editor"]["colnames"].text(5 + 70*(coords.mouseX + coords.Xshift), 10, tocompare),
                                {"clip": 70*(coords.mouseX + coords.Xshift) + ", 0, 68, 20"});
                            activeSquare();
                        string_command = "colnames(" + commobj["data_editor"]["dataset"] + ")[" + (coords.mouseX + coords.Xshift + 1) + "] <- " + ((tocompare == "")?"\"\"":("\"" + tocompare + "\""));
                        Rcommand.command = string_command; 
                        talkToR();
                    }
                    tobe.remove();
                    input.cover.toFront();
                    tasta = "enter";
                });
                colsrect.show();
                colsrect.toFront();
            });
            var rownamescover = papers["data_editor"]["rownames"].rect(0, 0, 70, 20*info["data"][dataset].nrows)
            .attr({fill: "#ffffff", stroke: "none", "fill-opacity": "0"})
            .click(function(event) {
                var coords = getCoords(event);
                bodyrect.hide();
                rowsrect.remove();
                colsrect.hide();
                papers["data_editor"]["topleft"].bodyrect_show = false;
                papers["data_editor"]["topleft"].rowsrect_show = true;
                papers["data_editor"]["topleft"].colsrect_show = false;
                rowsrect = papers["data_editor"]["rownames"].rect(1, 20*(coords.mouseY + coords.Yshift) + 1, 68, 18).attr({"stroke-width": 1.3});  
                papers["data_editor"]["topleft"].rowsrect = 20*(coords.mouseY + coords.Yshift) + 1;
                if ($("#data_editor_dataset_border").length) {
                    $("#data_editor_dataset_border").remove();
                }
            })
            .dblclick(function(event) {
                var coords = getCoords(event);
                rowsrect.hide();
                temp = info["data"][dataset].rownames[coords.mouseY + coords.Yshift];
                tobe = sat(papers["data_editor"]["rownames"].text(0, 0, temp));
                papers["data_editor"]["rownames"].inlineTextEditing(tobe);
                input = tobe.inlineTextEditing.startEditing(
                    0,
                    20*coords.mouseY - coords.scrollY + 20 + 20 - 1*(navigator.browserType == "Firefox"), 
                    70, 
                    20,
                    "whatever",
                    "#f2f2f2");
                input.cover = this;
                input.addEventListener("blur", function(e) {
                    tobe.inlineTextEditing.stopEditing(tasta);
                    tocompare = tobe.attr("text");
                    if (temp != tocompare) {
                        papers["data_editor"]["rownames"].rect(0, 20*(coords.mouseY + coords.Yshift), 70, 20)
                        .attr({fill: "#f2f2f2", stroke: "#d7d7d7"});
                        sat(papers["data_editor"]["rownames"].text(5, 10 + 20*(coords.mouseY + coords.Yshift), tocompare).toFront(),
                            {"clip": "0, " + 20*(coords.mouseY + coords.Yshift) + ", 68, 20"});
                        activeSquare();
                        string_command = "rownames(" + commobj["data_editor"]["dataset"] + ")[" + (coords.mouseY + coords.Yshift + 1) + "] <- " + ((tocompare == "")?"\"\"":("\"" + tocompare + "\""));
                        Rcommand.command = string_command; 
                        talkToR();
                    }
                    tobe.remove();
                    input.cover.toFront();
                    tasta = "enter";
                })
                rowsrect.show();
                rowsrect.toFront();
            })
            var datacover = papers["data_editor"]["body"].rect(0, 0, 70*info["data"][dataset].ncols, 20*info["data"][dataset].nrows)
            .attr({fill: "#aedaca", stroke: "none", "fill-opacity": 0})
            .click(function(event) {
                var coords = getCoords(event);
                bodyrect.remove();
                rowsrect.hide();
                colsrect.hide();
                papers["data_editor"]["topleft"].bodyrect_show = true;
                papers["data_editor"]["topleft"].rowsrect_show = false;
                papers["data_editor"]["topleft"].colsrect_show = false;
                bodyrect = papers["data_editor"]["body"].rect(70*(coords.mouseX + coords.Xshift) + 1, 20*(coords.mouseY + coords.Yshift) + 1, 68, 18)
                .attr({"stroke-width": 1.3});
                papers["data_editor"]["topleft"].bodyrect[0] = 70*(coords.mouseX + coords.Xshift) + 1;
                papers["data_editor"]["topleft"].bodyrect[1] = 20*(coords.mouseY + coords.Yshift) + 1;
                if ($("#data_editor_dataset_border").length) {
                    $("#data_editor_dataset_border").remove();
                }
            })  
            .dblclick(function(event) {
                var coords = getCoords(event);
                bodyrect.hide();
                temp = "" + info["data"][dataset].theData[coords.mX][coords.mY];
                temp = (temp == "null")?"":temp;
                tobe = sat(papers["data_editor"]["body"].text(0, 0, temp));
                papers["data_editor"]["body"].inlineTextEditing(tobe);
                input = tobe.inlineTextEditing.startEditing(
                    70*coords.mouseX - coords.scrollX + 70 + 1, 
                    20*coords.mouseY - coords.scrollY + 20 + 20 + 1 - 1*(navigator.browserType == "Firefox"), 
                    70, 
                    20 - 2);
                input.addEventListener("blur", function(e) {
                    tobe.inlineTextEditing.stopEditing(tasta);
                    tocompare = tobe.attr("text");
                    if (temp != tocompare) {
                        if (tocompare.length == 0) {
                            tocompare = "NA";
                        }
                        else {
                            if (!isNumeric(tocompare)) {
                                tocompare = "\"" + tocompare + "\"";
                            }
                        }
                        papers["data_editor"]["body"].rect(70*(coords.mouseX + coords.Xshift), 20*(coords.mouseY + coords.Yshift), 70, 20)
                        .attr({fill: "#ffffff", stroke: "none"});
                        sat(papers["data_editor"]["body"].text(5 + 70*(coords.mouseX + coords.Xshift), 10 + 20*(coords.mouseY + coords.Yshift), tocompare));
                        activeSquare();
                        string_command = dataset + "[\"" + info["data"][dataset].rownames[coords.mouseY + coords.Yshift] + "\", \"" + info["data"][dataset].colnames[coords.mouseX + coords.Xshift] + "\"] <- " + tocompare;
                        Rcommand.command = string_command; 
                        talkToR();
                    }
                    tobe.remove();
                    gridset.toFront();
                    bodyrect.show();
                    bodyrect.toFront();
                    datacover.toFront();
                    tasta = "enter";
                });
            });
        }
        else {
            var colnamescover = papers["data_editor"]["colnames"].rect(0, 0, 70*(visiblecols + 1), 20)
                .attr({fill: "#ffffff", stroke: "none", "fill-opacity": 0})
                .click(function(event) {
                    if ($("#data_editor_dataset_border").length) {
                        $("#data_editor_dataset_border").remove();
                    }
                });
            var rownamescover = papers["data_editor"]["rownames"].rect(0, 0, 70, 20*(visiblerows + 1))
                .attr({fill: "#ffffff", stroke: "none", "fill-opacity": "0"})
                .click(function(event) {
                    if ($("#data_editor_dataset_border").length) {
                        $("#data_editor_dataset_border").remove();
                    }
                });
            var datacover = papers["data_editor"]["body"].rect(0, 0, 70*(visiblecols + 1), 20*(visiblerows + 1))
                .attr({fill: "#aedaca", stroke: "none", "fill-opacity": 0})
                .click(function(event) {
                    if ($("#data_editor_dataset_border").length) {
                        $("#data_editor_dataset_border").remove();
                    }
                });
        }
        activeSquare();
    }
}
function refresh_dirs() {
    dirsfilescopy = "";
    if (dirfile.dirs != null) {
        dirsfilescopy += dirfile.dirs.toString();
    }
    if (dirfile.files != null) {
        dirsfilescopy += dirfile.files.toString();
    }
    dirfilevisit = true;
    print_dirs();
    updatecounter = 0;
    dirfilist.value = 1 - dirfilist.value;
    Shiny.onInputChange("dirfilist", dirfilist);
    printIfDirsFilesChange();
}
function draw_load(paper) {
    sat(paper.text(20, 20, "Package:"));
    sat(paper.text(190, 20, "Dataset:"));
    sat(paper.text(34, 230, "Help"));
    paper.rect(34 - 19, 230 - 13, 70, 25).attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
         .click(function() {
             if (commobj.load.pkgdata != "") {
                 Rcommand.command = "?" + commobj.load.pkgdata;
                 talkToR();
             }
             else if (commobj.load.package != "") {
                 Rcommand.command = "?" + commobj.load.package;
                 talkToR();
             }
         });
    sat(paper.text(355, 230, "Load"));
    paper.rect(355 - 19, 230 - 13, 70, 25).attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
         .click(function() {
             if (commobj.load.pkgdata != "" && !spinnerstarted) {
                 $("#load").remove();
                 Rcommand.command = string_command;
                 talkToR();
             }
         });
}
function draw_import(paper) {
    paper.clear();
    var stx = 13;
    var sty = 10;
    sat(paper.text(stx + 5, sty + 15, "Separator:"));
    var radios = paper.radio({
        x: stx + 11,
        y: sty + 40,
        whichChecked: 0,
        labels: ["comma", "space", "tab", "other, please specify:"],
    });
    radios.cover[0].click(function() {
        commobj.read_table.sep = ",";
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            Shiny.onInputChange("read_table", commobj.read_table);
            checkIfDataLoadedInR();
        }
    });
    radios.cover[1].click(function() {
        commobj.read_table.sep = " ";
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            Shiny.onInputChange("read_table", commobj.read_table);
            checkIfDataLoadedInR();
        }
    });
    radios.cover[2].click(function() {
        commobj.read_table.sep = "tab";
        Shiny.onInputChange("read_table", commobj.read_table);
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            checkIfDataLoadedInR();
        }
    });
    radios.cover[3].click(function() {
        commobj.read_table.sep = other.attr("text");
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            Shiny.onInputChange("read_table", commobj.read_table);
            checkIfDataLoadedInR();
        }
    });
    var other = sat(paper.text(stx + 170, sty + 116, ""),
                    {"clip": (stx + 165) + ", " + (sty + 106) + ", 35, 20"});
    var other_rect = paper.rect(stx + 165, sty + 106, 37, 20, 3).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"});
    paper.inlineTextEditing(other);
    var other_clicked = false;
    other_rect.click(function(e) {
        var me = this;
        e.stopPropagation();
        var temp = other.attr("text");
        ovBox = this.getBBox();
        input = other.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
        input.addEventListener("blur", function() {
            other.inlineTextEditing.stopEditing(tasta);
            if (other.attr("text") != temp) {
                commobj.read_table.sep = other.attr("text");
                radios.moveTo(3);
                if (dirfile.filename != "") {
                    console_command("import");
                    tempdatainfo.nrows = 0;
                    Shiny.onInputChange("read_table", commobj.read_table);
                    checkIfDataLoadedInR();
                }
            }
            me.toFront();
            tasta = "enter";
        }, true);
    });
    other_rect.mouseover(function() {
        if (other_clicked) {
            this.attr({'cursor':'pointer'});
        }
    });
    other_rect.mouseout(function() {
        this.attr({'cursor':''});
    });
    sat(paper.text(stx + 140, sty + 15, "Decimal:"));
    var decimal = paper.radio({
        x: stx + 150,
        y: sty + 40,
        whichChecked: 0,
        labels: ["dot", "comma"],
    });
    decimal.cover[0].click(function() {
        commobj.read_table.dec = ".";
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            Shiny.onInputChange("read_table", commobj.read_table);
            checkIfDataLoadedInR();
        }
    });
    decimal.cover[1].click(function() {
        commobj.read_table.dec = ",";
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            Shiny.onInputChange("read_table", commobj.read_table);
            checkIfDataLoadedInR();
        }
    });
    var header = paper.checkBox({
        x: stx + 5,
        y: sty + 142,
        isChecked: commobj.read_table.header,
        label: "Column names in the file header",
    });
    header.cover.click(function() {
        commobj.read_table.header = header.isChecked;
        if (dirfile.filename != "") {
            console_command("import");
            tempdatainfo.nrows = 0;
            Shiny.onInputChange("read_table", commobj.read_table);
            checkIfDataLoadedInR();
        }
    });
    var row_names = paper.rect(stx + 5, sty + 173, 70, 20, 3).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"});
    row_names.text = sat(paper.text(stx + 10, sty + 183, commobj.read_table.row_names),
                        {"clip": (stx + 5) + ", " + (sty + 173) + ", 68, 20"});
    sat(paper.text(stx + 80, sty + 175, "No./name of column"));
    sat(paper.text(stx + 80, sty + 191, "containing row names"));
    paper.inlineTextEditing(row_names.text);
    row_names.click(function(e) {
        var me = this;
        e.stopPropagation();
        var temp = row_names.text.attr("text");
        ovBox = me.getBBox();
        input = row_names.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
        input.addEventListener("blur", function(e) {
            row_names.text.inlineTextEditing.stopEditing(tasta);
            if (row_names.text.attr("text") != temp) {
                commobj.read_table.row_names = row_names.text.attr("text");
                Shiny.onInputChange("read_table", commobj.read_table);
                if (dirfile.filename != "") {
                    console_command("import");
                    tempdatainfo.nrows = 0;
                    checkIfDataLoadedInR();
                }
            }
            me.toFront();
            tasta = "enter";
        }, true);
    });
    sat(paper.text(stx + 5, sty + 218, "Preview column names:"));
    sat(paper.text(stx + 5, sty + 342, "Encoding:"));
    var fileEncoding = sat(
        paper.text(stx + 85, sty + 342, commobj["read_table"].fileEncoding),
        {"clip": (stx + 80) + ", " + (sty + 332) + ", 110, 20"}
    );
    paper.inlineTextEditing(fileEncoding);
    var enc_rect = paper.rect(stx + 80, sty + 332, 116, 20)
        .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
        .click(function(e) {
            e.stopPropagation();
            var me = this;
            var ovBox = this.getBBox();
            input = fileEncoding.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function(e) {
                fileEncoding.inlineTextEditing.stopEditing(tasta);
                commobj["read_table"].fileEncoding = fileEncoding.attr("text");
                console_command("import");
                me.toFront();
                tasta = "enter";
            }, true);
        });
    sat(paper.text(stx + 241, sty + 15, "Directory:"));
    paper.stdir_text = sat(
        paper.text(stx + 309, sty + 15, ""),
        {"clip": (stx + 306) + ", " + (sty + 5) + ", 308, 20"}
    );
    paper.inlineTextEditing(paper.stdir_text);
    var stdir_rect = paper.rect(stx + 304, sty + 5, 314, 20)
        .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
        .click(function(e) {
            e.stopPropagation();
            var temp = paper.stdir_text.attr("text");
            var BBox = this.getBBox();
            input = paper.stdir_text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                paper.stdir_text.inlineTextEditing.stopEditing(tasta);
                if (paper.stdir_text.attr("text") != temp) {
                    if (paper.stdir_text.attr("text") == "") {
                        paper.glow.hide();
                        dirfile_chosen[2] = "";
                    }
                    else {
                        dirfile_chosen[0] = "dir";
                        dirfile_chosen[1] = "__stdir__";
                        dirfile_chosen[2] = paper.stdir_text.attr("text");
                        pathcopy = dirfile.filepath;
                        Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                        printDirsWhenPathChanges();
                    }
                }
                tasta = "enter";
            }, true);
        });
    paper.glow = stdir_rect.glow({
        color: "#ff0000",
        width: 4
    });
    paper.glow.hide();
    sat(paper.text(stx + 560, sty + 378.5, "Import"));
    var import_open = paper.rect(stx + 542, sty + 366, 75, 25)
        .attr({"stroke-width": 1.25, fill: "#ffffff", "fill-opacity": 0})
        .click(function(e) {
            e.stopPropagation();
            if (dirfile.filename != "" && tempdatainfo.rownames != "foobar" && !spinnerstarted) {
                dirfilist.refresh = false;
                function littleWait() {
                    updatecounter += 1;
                    if (updatecounter < 101) { 
                        if (responseR) {
                            updatecounter = 0; 
                            printWhenOutputChanges(); 
                            tempdatainfo.nrows = 0;
                            $("#import").remove();
                            refresh_dirs(); 
                        }
                        else {
                            setTimeout(littleWait, 50);
                        }
                    }
                }
                if (commobj.importobj != commobj.read_table) {
                    commobj.importobj = copy(commobj.read_table);
                }
                responseR = false;
                updatecounter = 0;
                Shiny.onInputChange("import", commobj.importobj);
                littleWait();
            }
            else {
                if (tempdatainfo.rownames == "error!") {
                    history[(histindex < history.length)?history.length:histindex] = string_command;
                    histindex = history.length;
                    var header = strwrap(string_command, cpl.result, "  ");
                    $("#tempdiv").remove();
                    $("#result_main").append("<span style='color:#932192'>> </span><span style='color:blue'>" + header + "</span><br><br>");
                    $("#result_main").append("<span style='color:red'>Error: this is not a valid dataset.</span><br><br>");
                    createCommandPromptInRconsole();
                }
            }
        })
    var objname = paper.checkBox({
        x: stx + 1,
        y: sty + 373,
        isChecked: commobj.read_table.nameit,
        label: "Assign",
    });
    objname.cover.click(function() {
            console.log(1);
        if (this.isChecked) {
            objname.label[0].attr({"text": "Assign to:"});
            if (paper.objnametext.attr("text") == "") {
                commobj.read_table.customname = false;
                paper.objnametext.attr({"text": dirfile.filename});
            }
            commobj.read_table.objname = paper.objnametext.attr("text");
            commobj.read_table.nameit = true;
            objnameset.show();
        }
        else {
            objname.label[0].attr({"text": "Assign"});
            commobj.read_table.objname = "";
            commobj.read_table.nameit = false;
            objnameset.hide();
        }
        console_command("import");
    });
    var objnameset = paper.set();
    paper.objnametext = sat(paper.text(stx + 100, sty + 378, commobj.read_table.objname),
                    {"clip": (stx + 95) + ", " + (sty + 368) + ", 250, 20"});
    var objname_rect = paper.rect(stx + 95, sty + 368, 250, 20, 3).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"});
    paper.inlineTextEditing(paper.objnametext);
    objname_rect.click(function(e) {
        var me = this;
        e.stopPropagation();
        var temp = paper.objnametext.attr("text");
        ovBox = this.getBBox();
        input = paper.objnametext.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
        input.addEventListener("blur", function() {
            paper.objnametext.inlineTextEditing.stopEditing(tasta);
            if (paper.objnametext.attr("text") != temp) {
                commobj.read_table.objname = paper.objnametext.attr("text").replace(/[^A-Za-z0-9]/g, '');
                if (isNumeric(commobj.read_table.objname[0])) {
                    commobj.read_table.objname = "x" + commobj.read_table.objname;
                }
                paper.objnametext.attr({"text": commobj.read_table.objname});
                if (commobj.read_table.objname == "") {
                    commobj.read_table.customname = false;
                    paper.objnametext.attr({"text": dirfile.filename});
                    commobj.read_table.objname = dirfile.filename;
                }
                else {
                    commobj.read_table.customname = true;
                }
                console_command("import");
            }
            me.toFront();
            tasta = "enter";
        }, true);
    });
    objnameset.push(paper.objnametext, objname_rect);
    if (objname.isChecked) {
        objname.label[0].attr({"text": "Assign to:"});
        commobj.read_table.objname = paper.objnametext.attr("text");
        objnameset.show();
    }
    else {
        objnameset.hide();
    }
    console_command("import");
    refresh_dirs();
}
function draw_export(paper) {
    if ($("#export").length) {
        paper.clear();
        var stx = 13;
        var sty = 10;
        sat(paper.text(stx + 5, sty + 15, "Separator:"));
        var radios = paper.radio({
            x: stx + 11,
            y: sty + 40,
            whichChecked: 0,
            labels: ["comma", "space", "tab", "other, please specify:"],
        });
        radios.cover[0].click(function() {
            commobj["export"].sep = ",";
            console_command("export");
        });
        radios.cover[1].click(function() {
            commobj["export"].sep = " ";
            console_command("export");
        });
        radios.cover[2].click(function() {
            commobj["export"].sep = "tab";
            console_command("export");
        });
        var other = sat(paper.text(stx + 170, sty + 116, ""),
                        {"clip": (stx + 165) + "," + (sty + 111) + ",35,20"});
        var other_rect = paper.rect(stx + 165, sty + 106, 37, 20, 3).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"});
        paper.inlineTextEditing(other);
        var other_clicked = false;
        other_rect.click(function(e) {
            e.stopPropagation();
            var me = this;
            ovBox = this.getBBox();
            input = other.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function(e) {
                other.inlineTextEditing.stopEditing(tasta);
                commobj["export"].sep = other.attr("text");
                radios.moveTo(3);
                me.toFront();
                tasta = "enter";
            }, true);
        });
        other_rect.mouseover(function() {
            if (other_clicked) {
                this.attr({'cursor':'pointer'});
            }
        });
        other_rect.mouseout(function() {
            this.attr({'cursor':''});
        });
        var header = paper.checkBox({
            x: stx + 5,
            y: sty + 145,
            isChecked: commobj["export"].header,
            label: "Write column names",
        });
        header.cover.click(function() {
            commobj["export"].header = header.isChecked;
            console_command("export");
            if (header.isChecked) {
                caseidset.show();
            }
            else {
                caseidset.hide();
            }
        });
        var caseidset = paper.set();
        caseidset.push(sat(paper.text(stx + 30, sty + 175, "Cases ID:")));
        var caseid = sat(paper.text(stx + 107, sty + 175, commobj["export"].caseid),
                         {"clip": (stx + 107) + ", " + (sty + 165) + ", 98, 20"});
        var caseid_rect = paper.rect(stx + 102, sty + 165, 100, 20, 3).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"});
        paper.inlineTextEditing(caseid);
        caseid_rect.click(function(e) {
            e.stopPropagation();
            var me = this;
            ovBox = this.getBBox();
            input = caseid.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function(e) {
                caseid.inlineTextEditing.stopEditing(tasta);
                commobj["export"].caseid = caseid.attr("text");
                console_command("export");
                me.toFront();
                tasta = "enter";
            }, true);
        });
        caseidset.push(caseid, caseid_rect);
        sat(paper.text(stx + 5, sty + 198, "Dataset:"));
        sat(paper.text(stx + 5, sty + 317, "New file:"));
        paper.newname = sat(paper.text(stx + 74, sty + 323, commobj["export"].filename),
                            {"clip": (stx + 69) + ", " + (sty + 313) + ", 188, 20"});
        var newname_rect = paper.rect(stx + 69, sty + 313, 190, 20, 3).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"});
        paper.inlineTextEditing(paper.newname);
        newname_rect.click(function(e) {
            e.stopPropagation();
            var me = this;
            ovBox = this.getBBox();
            input = paper.newname.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function(e) {
                paper.newname.inlineTextEditing.stopEditing(tasta);
                commobj["export"].filename = paper.newname.attr("text");
                if (dirfile.files.indexOf(commobj["export"].filename) >= 0) {
                    paper.ovr.showIt();
                }
                else {
                    paper.ovr.hideIt();
                }
                console_command("export");
                me.toFront();
                tasta = "enter";
            }, true);
        });
        paper.ovr = paper.checkBox({
            x: stx + 270,
            y: sty + 317,
            isChecked: true,
            label: "Overwrite?",
        });
        paper.ovr.hideIt();
        if (dirfile.files.indexOf(commobj["export"].filename) >= 0) {
            paper.ovr.showIt();
        }
        paper.ovr.cover.click(function() {
            paper.newname.attr({"text": ""});
            commobj["export"].filename = "";
            paper.ovr.check();
            paper.ovr.hideIt();
            console_command("export");
        });
        sat(paper.text(stx + 217, sty + 15, "Directory:"));
        paper.stdir_text = sat(paper.text(stx + 286, sty + 15, ""),
                            {"clip": (stx + 283) + ", " + (sty + 5) + ", 286, 20"});
        paper.inlineTextEditing(paper.stdir_text);
        var stdir_rect = paper.rect(stx + 281, sty + 5, 292, 20)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(e) {
                e.stopPropagation();
                var temp = paper.stdir_text.attr("text");
                var BBox = this.getBBox();
                input = paper.stdir_text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    paper.stdir_text.inlineTextEditing.stopEditing(tasta);
                    if (paper.stdir_text.attr("text") != temp) {
                        if (paper.stdir_text.attr("text") == "") {
                            paper.glow.hide();
                            dirfile_chosen[2] = "";
                        }
                        else {
                            dirfile_chosen[0] = "dir";
                            dirfile_chosen[1] = "__stdir__";
                            dirfile_chosen[2] = paper.stdir_text.attr("text");
                            pathcopy = dirfile.filepath;
                            Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                            printDirsWhenPathChanges();
                        }
                    }
                    tasta = "enter";
                }, true);
            });
        paper.glow = stdir_rect.glow({
            color: "#ff0000",
            width: 4
        });
        paper.glow.hide();
        sat(paper.text(stx + 515, sty + 323.5, "Export"));
        paper.rect(stx + 497, sty + 311, 75, 25)
            .attr({"stroke-width": 1.25, fill: "#ffffff", "fill-opacity": 0})
            .click(function(e) {
                e.stopPropagation();
                if (!spinnerstarted) {
                    console_command("export");
                    Rcommand.command = string_command;
                    talkToR();
                    $("#export").remove();
                }
            })
        refresh_dirs();
    }
}
function draw_saveRplot(paper) {
    if ($("#saveRplot").length) {
        paper.clear();
        var stx = 13;
        var sty = 10;
        if (plotopen) {
            $("#saveRplot_preview").css({
                "background-image": "url('css/images/plot.svg?" + new Date().getTime() + "')", 
                "background-size": "100% 100%"
            });
        }
        var noplot = sat(paper.text(stx + 55, sty + 125, "No R plot window"));
        var filetypes = ["PNG", "BMP", "JPEG", "TIFF", "SVG", "PDF"];
        sat(paper.text(stx + 10, sty + 253, "Type:"));
        var radios = paper.radio({
            x: stx + 60,
            y: sty + 253,
            whichChecked: 0,
            labels: filetypes,
            vertspace: [0, 25, 50, 0, 25, 50],
            horspace: [0, 0, 0, 75, 75, 75],
        });
        radios.cover[0].click(function() {
            commobj.saveRplot.type = "png";
            if (commobj.saveRplot.filename !== "") {
                paper.filename.text.attr({"text": commobj.saveRplot.filename + "." + commobj.saveRplot.type})
            }
            checkOvr();
        });
        radios.cover[1].click(function() {
            commobj.saveRplot.type = "bmp";
            if (commobj.saveRplot.filename !== "") {
                paper.filename.text.attr({"text": commobj.saveRplot.filename + "." + commobj.saveRplot.type})
            }
            checkOvr();
        });
        radios.cover[2].click(function() {
            commobj.saveRplot.type = "jpeg";
            if (commobj.saveRplot.filename !== "") {
                paper.filename.text.attr({"text": commobj.saveRplot.filename + "." + commobj.saveRplot.type})
            }
            checkOvr();
        });
        radios.cover[3].click(function() {
            commobj.saveRplot.type = "tiff";
            if (commobj.saveRplot.filename !== "") {
                paper.filename.text.attr({"text": commobj.saveRplot.filename + "." + commobj.saveRplot.type})
            }
            checkOvr();
        });
        radios.cover[4].click(function() {
            commobj.saveRplot.type = "svg";
            if (commobj.saveRplot.filename !== "") {
                paper.filename.text.attr({"text": commobj.saveRplot.filename + "." + commobj.saveRplot.type})
            }
            checkOvr();
        });
        radios.cover[5].click(function() {
            commobj.saveRplot.type = "pdf";
            if (commobj.saveRplot.filename !== "") {
                paper.filename.text.attr({"text": commobj.saveRplot.filename + "." + commobj.saveRplot.type})
            }
            checkOvr();
        });
        var checkOvr = function() {
            if (dirfile.files.indexOf(commobj.saveRplot.filename + "." + commobj.saveRplot.type) >= 0) {
                paper.ovr.showIt();
            }
            else {
                paper.ovr.hideIt();
            }
        }
        sat(paper.text(stx, sty + 342, "File name:"));
        var splitext = function(text) {
            return(result = {
                filename: (text.length > 1)?(copy(text).splice(0, text.length - 1).join(".")):(text[0]),
                extension: (text.length > 1)?(text[text.length - 1]):("")
            })
        }
        var filename = (commobj.saveRplot.filename == "")?"":(commobj.saveRplot.filename + "." + commobj.saveRplot.type);
        paper.filename = textbox(paper, {x: stx + 74, y: sty + 343, width: 155, height: 20, text: filename});
        paper.inlineTextEditing(paper.filename.text);
        paper.filename.rect.click(function(e) {
            e.stopPropagation();
            var me = this;
            ovBox = this.getBBox();
            input = paper.filename.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function(e) {
                paper.filename.text.inlineTextEditing.stopEditing(tasta);
                var checkext = splitext(paper.filename.text.attr("text").split("."));
                var position = filetypes.indexOf(checkext.extension.toUpperCase());
                if (checkext.extension !== commobj.saveRplot.type) {
                    if (position >= 0) {
                        radios.moveTo(position);
                        commobj.saveRplot.filename = checkext.filename;
                        commobj.saveRplot.type = checkext.extension.toLowerCase();
                    }
                    else {
                        commobj.saveRplot.filename = paper.filename.text.attr("text");
                        paper.filename.text.attr({"text": paper.filename.text.attr("text") + "." + commobj.saveRplot.type});
                    }
                }
                else {
                    commobj.saveRplot.filename = checkext.filename;
                    commobj.saveRplot.type = checkext.extension.toLowerCase();
                }
                if (dirfile.files.indexOf(commobj.saveRplot.filename + "." + commobj.saveRplot.type) >= 0) {
                    paper.ovr.showIt();
                }
                else {
                    paper.ovr.hideIt();
                }
                me.toFront();
                tasta = "enter";
            }, true);
        });
        paper.ovr = paper.checkBox({
            x: stx + 236,
            y: sty + 338,
            isChecked: true,
            label: "Overwrite?",
        });
        paper.ovr.cover.click(function() {
            paper.filename.text.attr({"text": ""});
            commobj.saveRplot.filename = "";
            paper.ovr.check();
            paper.ovr.hideIt();
        });
        sat(paper.text(stx + 227, sty + 15, "Directory:"));
        paper.stdir_text = sat(paper.text(stx + 297, sty + 15, ""),
                            {"clip": (stx + 294) + ", " + (sty + 5) + ", 307, 20"});
        paper.inlineTextEditing(paper.stdir_text);
        var stdir_rect = paper.rect(stx + 292, sty + 5, 313, 20)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(e) {
                e.stopPropagation();
                var temp = paper.stdir_text.attr("text");
                var BBox = this.getBBox();
                input = paper.stdir_text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    paper.stdir_text.inlineTextEditing.stopEditing(tasta);
                    if (paper.stdir_text.attr("text") != temp) {
                        if (paper.stdir_text.attr("text") == "") {
                            paper.glow.hide();
                            dirfile_chosen[2] = "";
                        }
                        else {
                            dirfile_chosen[0] = "dir";
                            dirfile_chosen[1] = "__stdir__";
                            dirfile_chosen[2] = paper.stdir_text.attr("text");
                            pathcopy = dirfile.filepath;
                            Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                            printDirsWhenPathChanges();
                        }
                    }
                    tasta = "enter";
                }, true);
            });
        paper.glow = stdir_rect.glow({
            color: "#ff0000",
            width: 4
        });
        paper.glow.hide();
        sat(paper.text(stx + 549, sty + 343.5, "Save"));
        var saveRplot_rect = paper.rect(stx + 527, sty + 331, 75, 25)
            .attr({"stroke-width": 1.25, fill: "#ffffff", "fill-opacity": 0})
            .click(function(e) {
                e.stopPropagation();
                if (commobj.saveRplot.filename !== "" && plotopen) {
                    commobj.saveRplot.counter += 1;
                    Shiny.onInputChange("saveRplot", commobj.saveRplot);
                    $("#saveRplot").remove();
                } 
            })
        refresh_dirs();
    }
}
function refresh_cols(dialogs, exclude) {
    if (dialogs == "all") {
        dialogs = ["load", "import", "export", "minimize", "tt", "calibrate", "recode", "xyplot"];
        if (exclude !== void 0) {
            dialogs = copy(dialogs, exclude);
        }
    }
    else {
        dialogs = [dialogs];
    }
    var datasets = [], tts = [];
    if (info !== void 0) { 
        if (info["data"] !== null) { 
            datasets = getKeys(info["data"]);
        }
        if (info["tt"] !== null) {
            tts = getKeys(info["tt"]);
        }
    }
    for (var i = 0; i < dialogs.length; i++) {
        if (dialogs[i] == "load") {
            if (packages != void 0 && $("#load").length) {
                var cols = getKeys(packages);
                if (cols.length == 1) {
                    click_col("load", "package", cols[0], others = false);
                }
                print_cols("load", "package",
                               {
                                    "selection": "single",
                                    "cols": cols,
                                    "selectable": ["all"]
                               });
                if (commobj[dialogs[i]].package !== "") {
                    cols = getKeys(packages[commobj.load.package]);
                    print_cols("load", "pkgdata",
                               {
                                    "selection": "single",
                                    "cols": cols,
                                    "selectable": ["all"]
                               });
                }
            }
        }
        if (dialogs[i] == "import") {
            if ($("#import").length) {
                print_cols("import", "cols",
                           {
                                "selection": "none",
                                "cols": tempdatainfo.colnames,
                                "selectable": ["all"] 
                           });
                if (tempdatainfo.colnames.length == 1) {
                    papers["import"]["cols"].setSize(240, 20);
                }
            }
        }
        if (dialogs[i] == "export") {
            if ($("#export").length) {
                if (datasets.length) {
                    var cols = getKeys(info["data"]);
                    if (cols.length == 1) {
                        click_col(dialogs[i], "dataset", cols[0], others = false);
                    }
                    print_cols("export", "dataset",
                               {
                                    "selection": "single",
                                    "cols": cols,
                                    "selectable": ["all"]
                               });
                }
                else {
                    papers["export"]["dataset"].clear();
                    $(papers["export"]["dataset"].canvas).height(0);
                }
            }
        }
        if (dialogs[i] == "data_editor") {
            if ($("#data_editor").length && datasets.length) {
                print_cols("data_editor", "dataset",
                           {
                                "selection": "single",
                                "cols": datasets,
                                "selectable": ["all"]
                           });
            }
        }
        if (dialogs[i] == "findRows") {
            if ($("#findRows").length) {
                if (tts.length) {
                    if (tts.length == 1) {
                        click_col(dialogs[i], "tt", tts[0], others = false);
                        commobj.findRows.obj = tts[0];
                    }
                    print_cols("findRows", "tt",
                              {
                                   "selection": "single",
                                   "cols": tts,
                                   "selectable": ["all"]
                              });
                }   
                else {
                    papers["findRows"]["tt"].clear();
                    $(papers["findRows"]["tt"].canvas).height(0);
                }
            }
        }
        if (dialogs[i] == "minimize") {
            if ($("#minimize").length) {
                if (datasets.length || tts.length) {
                    var cols = getKeys(info[commobj[dialogs[i]].source]);
                    if (cols.length == 1) {
                        click_col(dialogs[i], "dataset", cols[0], others = false);
                    }
                    print_cols("minimize", "dataset",
                               {
                                    "selection": "single",
                                    "cols": cols,
                                    "selectable": ["all"]
                               });
                    if (commobj[dialogs[i]].dataset !== "" && info[commobj[dialogs[i]].source] !== null) {
                        cols = info[commobj[dialogs[i]].source][commobj[dialogs[i]].dataset].colnames;
                        print_cols("minimize", "outcome",
                               {
                                    "selection": (commobj[dialogs[i]].source == "data")?"multiple":"none",
                                    "cols": cols,
                                    "selectable": ["numerics"]
                               });
                        print_cols("minimize", "conditions",
                               {
                                    "selection": (commobj[dialogs[i]].source == "data")?"multiple":"none",
                                    "cols": cols,
                                    "selectable": ["numerics"]
                               });
                    }
                    else {
                        papers["minimize"]["outcome"].clear();
                        $(papers["minimize"]["outcome"].canvas).height(0);
                        papers["minimize"]["conditions"].clear();
                        $(papers["minimize"]["conditions"].canvas).height(0);
                    }
                }
                else {
                    papers["minimize"]["dataset"].clear();
                    $(papers["minimize"]["dataset"].canvas).height(0);
                    papers["minimize"]["outcome"].clear();
                    $(papers["minimize"]["outcome"].canvas).height(0);
                    papers["minimize"]["conditions"].clear();
                    $(papers["minimize"]["conditions"].canvas).height(0);
                }
            }
        }
        if (dialogs[i] == "tt") {
            if ($("#tt").length) {
                if (datasets.length) {
                    if (datasets.length == 1) {
                        click_col(dialogs[i], "dataset", datasets[0]);
                    }
                    print_cols("tt", "dataset",
                               {
                                    "selection": "single",
                                    "cols": datasets,
                                    "selectable": ["all"]
                               });
                    if (commobj[dialogs[i]].dataset !== "") {
                        print_cols("tt", "outcome",
                                   {
                                        "selection": "single",
                                        "cols": info["data"][commobj[dialogs[i]].dataset].colnames,
                                        "selectable": ["numerics"]
                                   });
                        print_cols("tt", "conditions",
                                   {
                                        "selection": "multiple",
                                        "cols": info["data"][commobj[dialogs[i]].dataset].colnames,
                                        "selectable": ["numerics"]
                                   });
                    }
                    else {
                        papers["tt"]["outcome"].clear();
                        $(papers["tt"]["outcome"].canvas).height(0);
                        papers["tt"]["conditions"].clear();
                        $(papers["tt"]["conditions"].canvas).height(0);
                    }
                }
                else {
                    papers["tt"]["dataset"].clear();
                    $(papers["tt"]["dataset"].canvas).height(0);
                    papers["tt"]["outcome"].clear();
                    $(papers["tt"]["outcome"].canvas).height(0);
                    papers["tt"]["conditions"].clear();
                    $(papers["tt"]["conditions"].canvas).height(0);
                }
            }
        }
        if (dialogs[i] == "calibrate") {
            if ($("#calibrate").length) {
                if (datasets.length) {
                    if (datasets.length == 1) {
                        click_col(dialogs[i], "dataset", datasets[0]);
                    }
                    print_cols("calibrate", "dataset",
                           {
                                "selection": "single",
                                "cols": datasets,
                                "selectable": ["all"]
                           });
                    if (commobj[dialogs[i]].dataset !== "") {
                        print_cols("calibrate", "x",
                           {
                                "selection": "single",
                                "cols": info["data"][commobj[dialogs[i]].dataset].colnames,
                                "selectable": ["numerics"]
                           });
                    }
                    else {
                        papers["calibrate"]["x"].clear();
                        $(papers["calibrate"]["x"].canvas).height(0);
                    }
                }
                else {
                    papers["calibrate"]["dataset"].clear();
                    $(papers["calibrate"]["dataset"].canvas).height(0);
                    papers["calibrate"]["x"].clear();
                    $(papers["calibrate"]["x"].canvas).height(0);
                }
            }
        }
        if (dialogs[i] == "recode") {
            if ($("#recode").length) {
                if (datasets.length) {
                    if (datasets.length == 1) {
                        click_col(dialogs[i], "dataset", datasets[0]);
                    }
                    print_cols("recode", "dataset",
                           {
                                "selection": "single",
                                "cols": datasets,
                                "selectable": ["all"]
                           });
                    if (commobj[dialogs[i]].dataset !== "") {
                        print_cols("recode", "x",
                               {
                                    "selection": "single",
                                    "cols": info["data"][commobj[dialogs[i]].dataset].colnames,
                                    "selectable": ["all"]
                               });
                        print_cols("recode", "rules",
                               {
                                    "selection": "multiple",
                                    "cols": makeRules(commobj[dialogs[i]].oldv, commobj[dialogs[i]].newv),
                                    "selectable": ["all"]
                               });
                    }
                    else {
                        papers["recode"]["x"].clear();
                        $(papers["recode"]["x"].canvas).height(0);
                    }
                }
                else {
                    papers["recode"]["dataset"].clear();
                    $(papers["recode"]["dataset"].canvas).height(0);
                    papers["recode"]["x"].clear();
                    $(papers["recode"]["x"].canvas).height(0);
                }
            }
        }
        if (dialogs[i] == "xyplot") {
            if ($("#xyplot").length) {
                if (datasets.length) {
                    if (datasets.length == 1) {
                        click_col(dialogs[i], "dataset", datasets[0]);
                    }
                    print_cols("xyplot", "dataset",
                               {
                                    "selection": "single",
                                    "cols": datasets,
                                    "selectable": ["all"]
                               });
                    if (commobj[dialogs[i]].dataset !== "") {
                        print_cols("xyplot", "x",
                                   {
                                        "selection": "single",
                                        "cols": info["data"][commobj["xyplot"].dataset].colnames,
                                        "selectable": ["numerics", "calibrated"]
                                   });
                        print_cols("xyplot", "y",
                                   {
                                        "dialogs": "xyplot",
                                        "identifier": "y",
                                        "selection": "single",
                                        "cols": info["data"][commobj["xyplot"].dataset].colnames,
                                        "selectable": ["numerics", "calibrated"]
                                   });
                    }
                    else {
                        papers["xyplot"]["x"].clear();
                        $(papers["xyplot"]["x"].canvas).height(0);
                        papers["xyplot"]["y"].clear();
                        $(papers["xyplot"]["y"].canvas).height(0);
                    }
                }
                else {
                    papers["xyplot"]["dataset"].clear();
                    $(papers["xyplot"]["dataset"].canvas).height(0);
                    papers["xyplot"]["x"].clear();
                    $(papers["xyplot"]["x"].canvas).height(0);
                    papers["xyplot"]["y"].clear();
                    $(papers["xyplot"]["y"].canvas).height(0);
                }
            }
        }
    }
}
function draw_calib(paper) {
if ($("#calibrate").length) {
    paper.clear();
    paper.crfuz = 1;
    var thlabelsfuz = ["e", "c", "i"];
    var thlabelscrp = ["th1", "th2", "th3", "th4", "th5", "th6"]
    var thlabels = new Array(6);
    var thcovers = new Array(6);
    var increasing = false;
    paper.thsetter_frame = paper.rect(160, 183.5, 325, 120).attr({stroke: "#d0d0d0"});
    sat(paper.text(18, 23, "Dataset:"));
    sat(paper.text(18, 153, "Choose condition:"));
    var stx = 153, sty = 32;
    var crfuz = paper.radio({
        x: stx + 15,
        y: sty,
        whichChecked: 1*(commobj.calibrate.type == "fuzzy"),
        labels: ["crisp", "fuzzy"],
    });
    crfuz.cover[0].click(function() {
        paper.crfuz = 0;
        commobj.calibrate.type = "crisp";
        var bell = sbell.whichChecked == 1;
        changeLabels();
        for (var i = 0; i < 6; i++) {
            if (bell) {
                commobj.calibrate.thscopyfuzb[i] = ths[i].attr("text");
            }
            else {
                commobj.calibrate.thscopyfuzs[i] = ths[i].attr("text");
            }
            commobj.calibrate.thresholds[i] = commobj.calibrate.thscopycrp[i];
            commobj.calibrate.thnames[i] = "th" + i;
            if (i >= commobj.calibrate.nth) {
                thsets[i].hide();
            }
            ths[i].attr({"text": commobj.calibrate.thresholds[i]});
        }
        showCrisp();
        console_command("calibrate");
        if (paper.findth.isChecked && commobj.calibrate.x != "") {
            thinfoToR();
        }
        else {
            drawPointsAndThresholds();
        }
    });
    crfuz.cover[1].click(function() {
        paper.crfuz = 1;
        commobj.calibrate.type = "fuzzy";
        var bell = sbell.whichChecked == 1;
        changeLabels();
        for (var i = 0; i < 6; i++) {
            commobj.calibrate.thscopycrp[i] = ths[i].attr("text");
            if (bell) {
                commobj.calibrate.thresholds[i] = commobj.calibrate.thscopyfuzb[i];
            }
            else {
                commobj.calibrate.thresholds[i] = commobj.calibrate.thscopyfuzs[i];
            }
            ths[i].attr({"text": commobj.calibrate.thresholds[i]});
        }
        for (var i = 0; i < ((sbell.whichChecked == 0)?3:6); i++) {
            commobj.calibrate.thnames[i] = thlabels.sub(i) + ((sbell.whichChecked == 0)?"":((i < 3)?1:2));
        }
        showFuzzy();
        console_command("calibrate");
    });
    var thno = paper.counter({"x": stx + 180, "y": sty + 65, "startval": commobj.calibrate.nth, "maxval": 6, "width": 18, "textheight": getTextHeight("0", 10)});
    thno.label({"x": -33, "y": 0, "label": "Number of thresholds:", "anchor": "end"});
    thno.update = function() {
        commobj.calibrate.nth = thno.value;
        commobj.calibrate.thresholds = ["", "", "", "", "", ""];
        for (var i = 0; i < commobj.calibrate.nth; i++) {
            commobj.calibrate.thresholds[i] = ths[i].attr("text");
        }
        showths();
        console_command("calibrate");
        if (paper.findth.isChecked) {
            thinfoToR();
        }
        else {
            drawPointsAndThresholds();
        }
    }
    thno.down.click(function() {
        thno.update();
    });
    thno.up.click(function() {
        thno.update();
    });
    function showths() {
        for (var i = 1; i < 6; i++) {
            if (i < commobj.calibrate.nth) {
                thsets[i].show();
            }
            else {
                thsets[i].hide();
            }
        }
    }
    paper.findth = paper.checkBox({
        x: stx + 109,
        y: sty - 5,
        isChecked: commobj.calibrate.findth,
        label: "find thresholds",
    });
    paper.findth.cover.click(function() {
        commobj.calibrate.findth = paper.findth.isChecked;
        if (getKeys(colclicks).indexOf("calibrate") >= 0) {
            commobj.calibrate.x = getTrueKeys(colclicks.calibrate.x)[0];
            if (commobj.calibrate.x === void 0) {
                commobj.calibrate.x = "";
            }
            if (commobj.calibrate.findth && commobj.calibrate.x != "") {
                thinfoToR();
            }
        }
    });
    var jitter = paper.checkBox({
        x: stx + 109, 
        y: sty + 20,
        isChecked: commobj.calibrate.jitter,
        label: "jitter points",
    });
    jitter.cover.click(function() {
        commobj.calibrate.jitter = jitter.isChecked;
        if (commobj.calibrate.x != "") {
            drawPointsAndThresholds();
        }
    });
    var logistic = paper.checkBox({
        x: stx + 9,
        y: sty + 58,
        isChecked: commobj.calibrate.logistic,
        label: "logistic",
    });
    logistic.cover.click(function() {
        commobj.calibrate.logistic = logistic.isChecked;
        if (logistic.isChecked) {
            sbell.moveTo(0);
            for (var i = 0; i < 3; i++) {
                commobj.calibrate.thresholds[i] = ths[i].attr("text");
                commobj.calibrate.thnames[i] = thlabels.sub(i);
            }
        }
        else {
            idm.hide();
        }
        showFuzzy();
        console_command("calibrate")
    });
    var idm = paper.set();
    idm.push(sat(paper.text(stx + 135, sty + 63, "degree of membership")));
    var idmtext = sat(paper.text(stx + 85, sty + 63, commobj.calibrate.idm));
    idm.push(idmtext);
    idm.push(sat(paper.rect(stx + 80, sty + 53, 50, 20, 3))
        .click(function(e) {
            e.stopPropagation();
            var me = this;
            var BBox = this.getBBox();
            input = idmtext.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                idmtext.inlineTextEditing.stopEditing(tasta);
                commobj.calibrate.idm = idmtext.attr("text");
                me.toFront();
                tasta = "enter";
                showFuzzy();
                console_command("calibrate");
            });
        }));
    paper.inlineTextEditing(idmtext);
    var ecdf = paper.checkBox({
        x: stx + 9,
        y: sty + 83,
        isChecked: commobj.calibrate.ecdf,
        label: "ecdf",
    });
    ecdf.cover.click(function() {
        commobj.calibrate.ecdf = ecdf.isChecked;
        if (commobj.calibrate.ecdf) {
            sbell.moveTo(0);
            thsets[3].hide();
            thsets[4].hide();
            thsets[5].hide();
            logistic.uncheck();
            commobj.calibrate.logistic = false;
            idm.hide();
            changeLabels();
            for (var i = 0; i < 3; i++) {
                commobj.calibrate.thresholds[i] = ths[i].attr("text");
                commobj.calibrate.thnames[i] = thlabels.sub(i);
            }
        }
        showFuzzy();
        console_command("calibrate");
    });
    var sbell = paper.radio({
        x: stx + 115,
        y: sty, 
        whichChecked: 1 - commobj.calibrate.end, 
        labels: ["s-shaped", "bell-shaped"],
    });
    sbell.cover[0].click(function() {
        showFuzzy();
        console_command("calibrate");
    });
    sbell.cover[1].click(function() {
        showFuzzy();
        console_command("calibrate");
    });
    var incdec = paper.radio({
        x: stx + 250, 
        y: sty, 
        whichChecked: 1 - commobj.calibrate.increasing, 
        labels: ["increasing", "decreasing"],
    });
    incdec.cover[0].click(function() {
        commobj.calibrate.increasing = true;
        changeLabels();
        for (var i = 0; i < 1*((sbell.whichChecked == 0)?3:6); i++) {
            commobj.calibrate.thresholds[i] = ths[i].attr("text");
            commobj.calibrate.thnames[i] = thlabels.sub(i) + ((sbell.whichChecked == 0)?"":((i < 3)?1:2));
        }
        showFuzzy();
        console_command("calibrate");
    });
    incdec.cover[1].click(function() {
        commobj.calibrate.increasing = false;
        changeLabels();
        for (var i = 0; i < 1*((sbell.whichChecked == 0)?3:6); i++) {
            commobj.calibrate.thresholds[i] = ths[i].attr("text");
            commobj.calibrate.thnames[i] = thlabels.sub(i) + ((sbell.whichChecked == 0)?"":((i < 3)?1:2))
        }
        showFuzzy();
        console_command("calibrate");
    });
    function changeLabels() {
        if (crfuz.whichChecked == 0) { 
            for (var i = 0; i < 6; i++) {
                thlabels[i].attr({"text": thlabelscrp[i]});
            }
        }
        else { 
            for (var i = 0; i < 6; i++) {
                if (incdec.whichChecked == 0) { 
                    thlabels[i].attr({"text": thlabelsfuz[(i < 3)?(i):(5 - i)] + ((sbell.whichChecked == 0)?"":((i < 3)?1:2))});
                }
                else { 
                    thlabels[i].attr({"text": thlabelsfuz[(i < 3)?(2 - i):(i - 3)] + ((sbell.whichChecked == 0)?"":((i < 3)?1:2))});
                }
                commobj.calibrate.thnames[i] = thlabels[i].attr("text");
            }
        }
    }
    var thsets = new Array(6);
    for (var i = 0; i < 6; i++) {
        thsets[i] = paper.set();
    }
    stx = 160;
    sty = 150;
    for (var i = 0; i < 6; i++) {
        thlabels[i] = sat(paper.text(stx + 5 + i*55, sty, thlabelscrp[i]));
        ths[i] = sat(paper.text(stx + 5 + i*55, sty + 20, (commobj.calibrate.type == "crisp")?commobj.calibrate.thscopycrp[i]:((sbell.whichChecked == 1)?commobj.calibrate.thscopyfuzb[i]:commobj.calibrate.thscopyfuzs[i])),
                     {"clip": (stx + i*55) + "," + (sty + 10) + ", 47, 20"});
        thcovers[i] = sat(paper.rect(stx + i*55, sty + 10, 50, 20, 3));
        thcovers[i].i = i;
        thcovers[i].click(function(e) {
            var me = this;
            e.stopPropagation();
            var temp = ths[me.i].attr("text");
            var tobe = sat(paper.text(-10, -10, temp));
            paper.inlineTextEditing(tobe);
            var BBox = me.getBBox();
            input = tobe.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.i = me.i;
            input.addEventListener("blur", function(e) {
                tobe.inlineTextEditing.stopEditing(tasta);
                if (temp != tobe.attr("text")) {
                    var finaltext = tobe.attr("text");
                    if ($.isNumeric(finaltext)) {
                        finaltext = 1*finaltext;
                    }
                    else {
                        finaltext = "";
                    }
                    ths[me.i].attr({"text": finaltext});
                    commobj.calibrate.thresholds[me.i] = finaltext;
                    var sort = false;
                    if (paper.crfuz == 0) {
                        for (var j = commobj.calibrate.nth; j < 6; j++) {
                            commobj.calibrate.thresholds[j] = "";
                        }
                        sort = true;
                    }
                    else {
                        var filled = 0;
                        for (var j = 0; j < 6; j++) {
                            if (commobj.calibrate.thresholds[j] != "") {
                                filled += 1;
                            }
                        }
                        sort = (sbell.whichChecked == 1)?(filled == 6):(filled == 3);
                    }
                    if (sort) {
                        commobj.calibrate.thresholds = sortArray(commobj.calibrate.thresholds);
                        for (var j = 0; j < 6; j++) {
                            ths[j].attr({"text": commobj.calibrate.thresholds[j]});
                        }
                    }
                    if (paper.crfuz == 0) { 
                        commobj.calibrate.thscopycrp = copy(commobj.calibrate.thresholds);
                        paper.findth.uncheck();
                        commobj.calibrate.findth = false;
                        drawPointsAndThresholds();
                    }
                    else { 
                        if (sbell.whichChecked == 1) {
                            commobj.calibrate.thnames[me.i] = thlabels.sub(me.i) + ((me.i < 3)?1:2);
                            commobj.calibrate.thscopyfuzb = copy(commobj.calibrate.thresholds);
                        }
                        else {
                            commobj.calibrate.thnames[me.i] = thlabels.sub(me.i);
                            commobj.calibrate.thscopyfuzs = copy(commobj.calibrate.thresholds);
                        }
                        showFuzzy(); 
                    }
                    console_command("calibrate");
                }
                me.toFront();
                tasta = "enter";
                tobe.remove();
            });
        });
        thsets[i].push(ths[i], thlabels[i], thcovers[i]);
    }
    thlabels[0].attr({"text": thlabelscrp[0]});
    thlabels.sub = function(x) {
        return(thlabels[x].attr("text").substring(0, 1))
    }
    var abset = paper.set();
    abset.push(sat(paper.text(stx + 215, sty - 43, "Shape form:"), {"anchor": "end"}));
    abset.push(sat(paper.text(stx + 225, sty - 55, "above")));
    abset.push(sat(paper.text(stx + 225, sty - 30, "below")));
    var avalue = sat(paper.text(stx + 279.5, sty - 55, "1"), {"clip": (stx + 275) + "," + (sty - 65) + ", 47, 20"});
    var bvalue = sat(paper.text(stx + 279.5, sty - 30, "1"), {"anchor": "start"});
    abset.push(bvalue, avalue);
    abset.push(sat(paper.rect(stx + 275, sty - 65, 50, 20, 3))
        .click(function(e) {
            e.stopPropagation();
            var me = this;
            var temp = avalue.attr("text");
            var tobe = sat(paper.text(-10, -10, temp));
            paper.inlineTextEditing(tobe);
            var BBox = me.getBBox();
            input = tobe.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                tobe.inlineTextEditing.stopEditing(tasta);
                var finaltext = tobe.attr("text");
                if ($.isNumeric(finaltext)) {
                    finaltext = 1*finaltext;
                    if (finaltext <= 0) {
                        finaltext = temp;
                    }
                }
                else {
                    finaltext = temp;
                }
                avalue.attr({"text": finaltext});
                commobj.calibrate.above = finaltext;
                me.toFront();
                tasta = "enter";
                showFuzzy(); 
                console_command("calibrate");
            });
        }));
    abset.push(sat(paper.rect(stx + 275, sty - 40, 50, 20, 3))
        .click(function(e) {
            e.stopPropagation();
            var me = this;
            var temp = bvalue.attr("text");
            var tobe = sat(paper.text(-10, -10, temp));
            paper.inlineTextEditing(tobe);
            var BBox = this.getBBox();
            input = tobe.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                tobe.inlineTextEditing.stopEditing(tasta);
                var finaltext = tobe.attr("text");
                if ($.isNumeric(finaltext)) {
                    finaltext = 1*finaltext;
                    if (finaltext <= 0) {
                        finaltext = temp;
                    }
                }
                else {
                    finaltext = temp;
                }
                bvalue.attr({"text": finaltext});
                commobj.calibrate.below = finaltext;
                me.toFront();
                tasta = "enter";
                showFuzzy(); 
                console_command("calibrate");
            });
        }));
    function showCrisp() {
        thno.showIt();
        paper.thsetter_frame.show();
        thsetter_content.show();
        paper.findth.showIt();
        logistic.hideIt();
        ecdf.hideIt();
        incdec.hideIt();
        sbell.hideIt();
        idm.hide();
        abset.hide();
        jitter.showIt();
        showths();
        changeLabels();
        drawPointsAndThresholds();
    }
    function showFuzzy() {
        var bell = sbell.whichChecked == 1;
        commobj.calibrate.end = !bell;
        logistic.showIt();
        ecdf.showIt();
        if (bell) {
            thsets[3].show();
            thsets[4].show();
            thsets[5].show();
            commobj.calibrate.logistic = false;
            commobj.calibrate.ecdf = false;
            logistic.uncheck();
            ecdf.uncheck();
            idm.hide();
            abset.show();
        }
        else {
            thsets[3].hide();
            thsets[4].hide();
            thsets[5].hide();
            if (logistic.isChecked) {
                idm.show();
                abset.hide();
                ecdf.uncheck();
            }
            else {
                idm.hide();
                if (ecdf.isChecked) {
                    abset.hide();
                }
                else {
                    abset.show();
                }
            }
        }
        changeLabels();
        for (var i = 0; i < 6; i++) {
            commobj.calibrate.thresholds[i] = (crfuz.whichChecked == 0)?commobj.calibrate.thscopycrp[i]:((sbell.whichChecked == 0)?commobj.calibrate.thscopyfuzs[i]:commobj.calibrate.thscopyfuzb[i]);
            ths[i].attr({"text": commobj.calibrate.thresholds[i]});
        }
        incdec.showIt();
        sbell.showIt();
        paper.findth.hideIt();
        thno.hideIt();
        jitter.hideIt();
        thsets[1].show();
        thsets[2].show();
        if (colclicks.calibrate !== void 0) {
            commobj.calibrate.x = getTrueKeys(colclicks.calibrate.x)[0];
            if (commobj.calibrate.x === void 0) {
                commobj.calibrate.x = "";
            }
        }
        if (commobj.calibrate.x != "") {
            if (commobj.calibrate.type == "fuzzy") {
                if (checkfuzzy()) {
                    commobj.calibrate.counter += 1;
                    commobj.calibrate.thsetter = true;
                    Shiny.onInputChange("calibrate", commobj.calibrate);
                }
                else {
                    poinths.fuzzyvals = new Array();
                    drawPointsAndThresholds();
                }
            }
            else {
                drawPointsAndThresholds(); 
            }
        }
        else {
            drawPointsAndThresholds(); 
        }
    }
    if (commobj.calibrate.type == "crisp") {
        showCrisp();
    }
    else {
        showFuzzy();
    }
    changeLabels();
    stx = 13, sty = 265;
    var newcond = paper.checkBox({
        x: stx + 3,
        y: sty + 75,
        isChecked: !commobj.calibrate.same,
        label: "into new condition",
    });
    newcond.label[0].remove();
    newcond.label = new Array(2);
    newcond.label[0] = sat(paper.text(stx + 27, sty + 73, "calibrate into"));
    newcond.label[1] = sat(paper.text(stx + 27, sty + 88, "new condition"));
    newcond.cover.click(function() {
        commobj.calibrate.same = !newcond.isChecked;
        if (newcond.isChecked) {
            newname.show();
        }
        else {
            newname.hide();
        }
        console_command("calibrate");
    });
    var newname = paper.set();
    var newnametext = sat(paper.text(stx + 130, sty + 80, commobj.calibrate.newvar),
                          {"clip": (stx + 125) + "," + (sty + 70) + ", 147, 20"});
    newname.push(newnametext);
    newname.push(sat(paper.rect(stx + 125, sty + 70, 150, 20, 3))
        .click(function(e) {
            var me = this;
            e.stopPropagation();
            var BBox = this.getBBox();
            input = newnametext.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                newnametext.inlineTextEditing.stopEditing(tasta);
                commobj.calibrate.newvar = newnametext.attr("text");
                if (isNumeric(commobj.calibrate.newvar[0])) {
                    commobj.calibrate.newvar = "x" + commobj.calibrate.newvar;
                }
                newnametext.attr({"text": commobj.calibrate.newvar});
                me.toFront();
                tasta = "enter";
                console_command("calibrate");
            });
        }));
    paper.inlineTextEditing(newnametext);
    if (commobj.calibrate.same) {
        newname.hide();
    }
    paper.Run = paper.set();
    paper.Run.push(sat(paper.text(423, sty + 80, "Run")));
    paper.Run.push(paper.rect(421 - 20, sty + 67, 70, 25)
    .attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
    .click(function() {
        focus = "calibrate";
        if (info["data"][commobj["calibrate"].dataset].rownames != "" && !spinnerstarted) {
            commobj.calibrate.x = getTrueKeys(colclicks.calibrate.x)[0];
            if (commobj.calibrate.x === void 0) {
                commobj.calibrate.x = "";
            }
            if (commobj.calibrate.x != "" && checkfuzzy()) {
                Rcommand.command = string_command;
                talkToR();
            }
        }
    }));
    if (getKeys(colclicks).indexOf("calibrate") >= 0) {
        if (commobj.calibrate.findth && commobj.calibrate.x != "") {
            drawPointsAndThresholds();
        }
    }
} 
} 
function checkfuzzy() {
    var fuzcheck = true;
    if (commobj.calibrate.type == "fuzzy") {
        if (commobj.calibrate.logistic) {
            if (commobj.calibrate.idm == "") {
                fuzcheck = false;
            }
            else {
                fuzcheck = fuzcheck && !isNaN(commobj.calibrate.idm);
            }
        }
        if (commobj.calibrate.below == "" || commobj.calibrate.above == "") {
            fuzcheck = false;
        }
        else {
            fuzcheck = fuzcheck && !isNaN(commobj.calibrate.below) && !isNaN(commobj.calibrate.above)
        }
        for (var i = 0; i < 1*((commobj.calibrate.end)?3:6); i++) {
            if (commobj.calibrate.thresholds[i] == "") {
                fuzcheck = false;
            }
        }
    }
    return(fuzcheck)
}
function draw_recode(paper) {
if ($("#recode").length) {
    paper.clear();
    paper.rules = new Array();
    paper.rules["oldv"] = "";
    paper.rules["newv"] = "";
    var stx = 13, sty = 10;
    sat(paper.text(stx + 5, sty + 13, "Dataset:"), {"text": 0});
    sat(paper.text(stx + 5, sty + 123, "Choose condition:"), {"text": 0});
    paper.newcond = paper.checkBox({
        x: stx + 1,
        y: sty + 280,
        isChecked: !commobj.recode.same,
        label: "into new condition",
    });
    paper.newcond.label[0].remove();
    paper.newcond.label = new Array(2);
    paper.newcond.label[0] = sat(paper.text(stx + 25, sty + 278, "recode into"), {"text": 0});
    paper.newcond.label[1] = sat(paper.text(stx + 25, sty + 293, "new condition"), {"text": 0});
    paper.newcond.cover.click(function() {
        commobj.recode.same = !paper.newcond.isChecked;
        if (paper.newcond.isChecked) {
            paper.newnameset.show();
        }
        else {
            paper.newnameset.hide();
        }
        console_command("recode");
    });
    paper.newnameset = paper.set();
    paper.newname_TB = textbox(paper, {x: stx + 128, y: sty + 285, width: 150, height: 20, text: commobj.recode.newvar});
    paper.inlineTextEditing(paper.newname_TB.text);
    paper.newname_TB.rect.click(function(e) {
        var me = this;
        e.stopPropagation();
        var temp = paper.newname_TB.text.attr("text");
        ovBox = this.getBBox();
        input = paper.newname_TB.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
        input.addEventListener("blur", function() {
            paper.newname_TB.text.inlineTextEditing.stopEditing(tasta);
            commobj.recode.newvar = paper.newname_TB.text.attr("text");
            if (isNumeric(commobj.recode.newvar[0])) {
                commobj.recode.newvar = "x" + commobj.recode.newvar;
            }
            paper.newname_TB.text.attr({"text": commobj.recode.newvar});
            me.toFront();
            tasta = "enter";
            console_command("recode");
        }, true);
    });
    paper.newnameset.push(paper.newname_TB.text, paper.newname_TB.rect);
    if (commobj.recode.same) {
        paper.newnameset.hide();
    }
    var stx = 163, sty = 24, vertspace = 38;
    paper.oldv = new Array();
    paper.newv = new Array();
    paper.oldv.texts = new Array()
    paper.oldv.covers = new Array();
    paper.newv.texts = new Array()
    paper.newv.covers = new Array();
    paper.text(stx - 6, sty, "Old value(s):").attr({"text-anchor": "start", "font-size": "14px"});
    paper.text(stx + 175, sty, "New value:").attr({"text-anchor": "start", "font-size": "14px"});
    paper.oldradio = paper.radio({
        x: stx,
        y: sty + vertspace, 
        whichChecked: -1, 
        labels: [
            "value",
            "lowest to",
            "\u00a0 \u00a0 \u00a0 \u00a0 \u00a0 \u00a0 to",
            "\u00a0 \u00a0 \u00a0 \u00a0 \u00a0 \u00a0 to highest",
            "missing",
            "all other values"
            ],
        vertspace: vertspace,
    });
    paper.oldradio.cover[0].click(function() {
        paper.rules.oldv = paper.oldv.texts.VALUE.attr("text");
    });
    paper.oldradio.cover[1].click(function() {
        paper.rules.oldv = paper.oldv.texts.LOWESTTO.attr("text");
    });
    paper.oldradio.cover[2].click(function() {
        paper.rules.oldv[0] = paper.oldv.texts.range.FROM.attr("text");
        paper.rules.oldv[1] = paper.oldv.texts.range.TO.attr("text");
    });
    paper.oldradio.cover[3].click(function() {
        paper.rules.oldv = paper.oldv.texts.TOHIGHEST.attr("text");
    });
    paper.oldradio.cover[4].click(function() {
        paper.rules.oldv = "missing";
    });
    paper.oldradio.cover[5].click(function() {
        paper.rules.oldv = "else";
    });
    paper.oldv.texts.VALUE  = paper.text(stx + 60, sty + vertspace, "");
    paper.oldv.covers.VALUE = paper.rect(stx + 55, sty + vertspace - 10, 40, 20, 3)
    .click(function(e) {
        e.stopPropagation();
        var me = this;
        paper.oldradio.moveTo(0);
        var BBox = this.getBBox();
        input = paper.oldv.texts.VALUE.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
        input.addEventListener("blur", function(e) {
            paper.oldv.texts.VALUE.inlineTextEditing.stopEditing(tasta);
            paper.rules.oldv = paper.oldv.texts.VALUE.attr("text");
            me.toFront();
            tasta = "enter";
        });
    });
    paper.inlineTextEditing(paper.oldv.texts.VALUE);
    paper.oldv.texts.LOWESTTO = paper.text(stx + 85, sty + 2*vertspace, "");
    paper.oldv.covers.LOWESTTO = paper.rect(stx + 80, sty + 2*vertspace - 10, 40, 20, 3)
    .click(function(e) {
        e.stopPropagation();
        var me = this;
        paper.oldradio.moveTo(1);
        var BBox = this.getBBox();
        input = paper.oldv.texts.LOWESTTO.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
        input.addEventListener("blur", function(e) {
            paper.oldv.texts.LOWESTTO.inlineTextEditing.stopEditing(tasta);
            paper.rules.oldv = paper.oldv.texts.LOWESTTO.attr("text");
            me.toFront();
            tasta = "enter";
        });
    });
    paper.inlineTextEditing(paper.oldv.texts.LOWESTTO);
    paper.oldv.texts.range = new Array();
    paper.oldv.covers.range = new Array();
    paper.oldv.texts.range.FROM  = paper.text(stx + 19, sty + 3*vertspace, "");
    paper.oldv.covers.range.FROM = paper.rect(stx + 14, sty + 3*vertspace - 10, 40, 20, 3)
    .click(function(e) {
        e.stopPropagation();
        paper.oldradio.moveTo(2);
        paper.rules.oldv = new Array(2);
        paper.rules.oldv[0] = paper.oldv.texts.range.FROM.attr("text");
        paper.rules.oldv[1] = paper.oldv.texts.range.TO.attr("text");
        var BBox = this.getBBox();
        input = paper.oldv.texts.range.FROM.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
        input.addEventListener("blur", function(e) {
            paper.oldv.texts.range.FROM.inlineTextEditing.stopEditing(tasta);
            paper.rules.oldv[0] = paper.oldv.texts.range.FROM.attr("text");
            tasta = "enter";
        });
    });
    paper.inlineTextEditing(paper.oldv.texts.range.FROM);
    paper.oldv.texts.range.TO = paper.text(stx + 85, sty + 3*vertspace, "");
    paper.oldv.covers.range.TO = paper.rect(stx + 80, sty + 3*vertspace - 10, 40, 20, 3)
    .click(function(e) {
        e.stopPropagation();
        var me = this;
        paper.oldradio.moveTo(2);
        paper.rules.oldv = new Array(2);
        paper.rules.oldv[0] = paper.oldv.texts.range.FROM.attr("text");
        paper.rules.oldv[1] = paper.oldv.texts.range.TO.attr("text");
        var BBox = this.getBBox();
        input = paper.oldv.texts.range.TO.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
        input.addEventListener("blur", function(e) {
            paper.oldv.texts.range.TO.inlineTextEditing.stopEditing(tasta);
            paper.rules.oldv[1] = paper.oldv.texts.range.TO.attr("text");
            me.toFront();
            tasta = "enter";
        });
    });
    paper.inlineTextEditing(paper.oldv.texts.range.TO);
    paper.oldv.texts.TOHIGHEST = paper.text(stx + 19, sty + 4*vertspace, "");
    paper.oldv.covers.TOHIGHEST = paper.rect(stx + 14, sty + 4*vertspace - 10, 40, 20, 3)
    .click(function(e) {
        e.stopPropagation();
        var me = this;
        paper.oldradio.moveTo(3);
        var BBox = this.getBBox();
        input = paper.oldv.texts.TOHIGHEST.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
        input.addEventListener("blur", function(e) {
            paper.oldv.texts.TOHIGHEST.inlineTextEditing.stopEditing(tasta);
            paper.rules.oldv = paper.oldv.texts.TOHIGHEST.attr("text");
            me.toFront();
            tasta = "enter";
        });
    });
    paper.inlineTextEditing(paper.oldv.texts.TOHIGHEST);
    sat(paper.oldv.texts, {"clip": paper.oldv.covers});
    sat(paper.oldv.covers);
    paper.path([ 
        ["M", stx + 140, sty + 15],
        ["L", stx + 140, sty + 240]
    ]).attr({stroke: "#a0a0a0"});
    paper.newradio = paper.radio({
        x: stx + 180, 
        y: sty + vertspace,
        whichChecked: -1,
        labels: ["value", "missing", "copy old value(s)"],
        vertspace: vertspace - 10,
    });
    paper.newradio.cover[0].click(function() {
        paper.rules.newv = paper.newv.texts.VALUE.attr("text");;
    });
    paper.newradio.cover[1].click(function() {
        paper.rules.newv = "missing";
    });
    paper.newradio.cover[2].click(function() {
        paper.rules.newv = "copy";
    });
    paper.newv.texts.VALUE  = paper.text(stx + 61 + 180, sty + vertspace, "");
    paper.newv.covers.VALUE = paper.rect(stx + 56 + 180, sty + vertspace - 10, 40, 20, 3)
    .attr({stroke: '#a0a0a0', 'stroke-width': 1, fill: "#ffffff", "fill-opacity": 0})
    .click(function(e) {
        e.stopPropagation();
        var me = this;
        paper.newradio.moveTo(0);
        var BBox = this.getBBox();
        input = paper.newv.texts.VALUE.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
        input.addEventListener("blur", function(e) {
            paper.newv.texts.VALUE.inlineTextEditing.stopEditing(tasta);
            paper.rules.newv = paper.newv.texts.VALUE.attr("text");
            me.toFront();
            tasta = "enter";
        });
    });
    paper.inlineTextEditing(paper.newv.texts.VALUE);
    sat(paper.newv.texts, {"clip": paper.newv.covers});
    sat(paper.newv.covers);
    paper.add = sat(paper.text(stx + 188, sty + 142, "Add"), {"size": 12, "anchor": "middle"});
    sat(paper.text(stx + 225, sty + 142, "Remove"), {size: 12});
    sat(paper.text(stx + 293, sty + 142, "Clear"), {size: 12});
    sat(paper.rect(stx + 163, sty + 132, 50, 20), {"sw": 1}) 
    .click(function() {
        var rule;
        if (all(paper.rules.oldv, " != \"\"") && paper.rules.newv != "") {
            if (paper.oldradio.whichChecked == 1) { 
                rule = "lo:" + paper.rules.oldv;
            }
            else if (paper.oldradio.whichChecked == 2) {
                rule = paper.rules.oldv[0] + ":" + paper.rules.oldv[1];
            }
            else if (paper.oldradio.whichChecked == 3) {
                rule = paper.rules.oldv + ":hi";
            }
            else {
                rule = paper.rules.oldv;
            }
            if (paper.add.attr("text") == "Change") {
                var selected = getTrueKeys(colclicks.recode.rules);
                if (selected.length == 1) {
                    for (var i = 0, len = commobj.recode.oldv.length; i < len; i++) {
                        if (commobj.recode.oldv[i] + "=" + commobj.recode.newv[i] == selected) {
                            commobj.recode.oldv[i] = rule;
                            commobj.recode.newv[i] = paper.rules.newv;
                        }
                    }
                }
            }
            else {
                commobj.recode.oldv.push(rule);
                commobj.recode.newv.push(paper.rules.newv);
            }
            if (colclicks.recode.rules != void 0) {
                unselect(colclicks, "recode", "rules");
            }
            eraseRecodeValues(paper);
            print_cols("recode", "rules",
                       {
                            "dialog": "recode",
                            "identifier": "rules",
                            "selection": "multiple",
                            "cols": makeRules(commobj.recode.oldv, commobj.recode.newv),
                            "selectable": ["all"]
                       });
        }
        console_command("recode");
        paper.add.attr({"text": "Add"});
    });
    sat(paper.rect(stx + 218, sty + 132, 60, 20), {"sw": 1}) 
    .click(function() {
        paper.add.attr({"text": "Add"});
        var selected = getTrueKeys(colclicks.recode.rules);
        for (var i = 0; i < selected.length; i++) {
            var idx = commobj.recode.oldv.indexOf(selected[i].split("=")[0]);
            commobj.recode.oldv.splice(idx, 1);
            commobj.recode.newv.splice(idx, 1);
            deleteRule(colclicks, selected[i]);
        }
        papers["recode"]["rules"].clear();
        eraseRecodeValues(paper);
        if (commobj.recode.oldv.length > 0) {
            print_cols("recode", "rules",
                       {
                            "dialog": "recode",
                            "identifier": "rules",
                            "selection": "multiple",
                            "cols": makeRules(commobj.recode.oldv, commobj.recode.newv),
                            "selectable": ["all"]
                       });
        }
        console_command("recode");
    });
    sat(paper.rect(stx + 283, sty + 132, 50, 20), {"sw": 1}) 
    .click(function() {
        paper.add.attr({"text": "Add"});
        papers["recode"]["rules"].clear();
        eraseRecodeValues(paper);
        commobj.recode.oldv = new Array();
        commobj.recode.newv = new Array();
        console_command("recode");
    });
    sat(paper.text(445, 295, "Run"));
    paper.rect(445 - 20, 295 - 13, 70, 25).attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
         .click(function() {
             if (info["data"][commobj["recode"].dataset].rownames != "" && !spinnerstarted) {
                 commobj.recode.x = getTrueKeys(colclicks.recode.x)[0];
                 if (commobj.recode.x === void 0) {
                     commobj.recode.x = "";
                 }
                 if (commobj.recode.x != "" && commobj.recode.newv.length > 0) {
                     Rcommand.command = string_command;
                     talkToR();
                 }
             }
         });
} 
} 
function draw_xyplot(paper) {
if ($("#xyplot").length) {
    var scale = 1;
    if (paper.scale != void 0) {
        scale = paper.scale;
    }
    if (paper.th != void 0) {
        paper.th.remove();
    }
    var labelRotation = 0;
    if (paper.labelRotation != void 0) {
        labelRotation = paper.labelRotation;
    }
    var randomjitter = new Array();
    if (paper.randomjitter != void 0) {
        if (getKeys(paper.randomjitter).length > 0) {
            randomjitter = paper.randomjitter;
        }
    }
    paper.clear();
    paper.xyplotdata = copy(xyplotdata);
    paper.randomjitter = randomjitter;
    paper.labelRotation = labelRotation;
    var stx = 13, sty = 10;
    sat(paper.text(stx, sty + 10, "Dataset:"));
    sat(paper.text(stx, sty + 125, "Condition X:"));
    sat(paper.text(stx, sty + 240, "Outcome Y:"));
    paper.scale = scale;
    paper.sx = 230;
    paper.sy = 20;
    paper.dim = 480; 
    paper.offset = 8;
    paper.rdim = paper.dim - 2*paper.offset;
    paper.negx = paper.checkBox({
        x: stx + 89,
        y: sty + 120,
        isChecked: commobj.xyplot.negx,
        label: "negate",
    });
    paper.negy = paper.checkBox({
        x: stx + 89,
        y: sty + 235,
        isChecked: commobj.xyplot.negy,
        label: "negate",
    });
    paper.index = 0;
    var powersof2 = 2;
    paper.negy.cover.click(function() {
        commobj.xyplot.negy = this.isChecked;
        powersof2 = Math.pow(2, (this.isChecked)?3:0) + Math.pow(2, (paper.negx.isChecked)?2:0);
        paper.index = [2, 5, 9, 12].indexOf(powersof2);
        paper.x = commobj.xyplot.x;
        paper.y = commobj.xyplot.y;
        scaleplot(paper);
        createLabels(paper);
    });
    paper.negx.cover.click(function() {
        commobj.xyplot.negx = this.isChecked;
        powersof2 = Math.pow(2, (paper.negy.isChecked)?3:0) + Math.pow(2, (this.isChecked)?2:0);
        paper.index = [2, 5, 9, 12].indexOf(powersof2);
        paper.x = commobj.xyplot.x;
        paper.y = commobj.xyplot.y;
        scaleplot(paper);
        createLabels(paper);
    });
    paper.sufnec = paper.radio({
        x: stx + 7,
        y: sty + 355,
        whichChecked: ["sufficiency", "necessity"].indexOf(commobj.xyplot.sufnec),
        labels: ["sufficiency", "necessity"],
    });
    paper.pof = paper.checkBox({
        x: stx + 1,
        y: sty + 405,
        isChecked: commobj.xyplot.pof,
        label: "parameters of fit",
    });
    paper.mdguides = paper.checkBox({
        x: stx + 1,
        y: sty + 430,
        isChecked: commobj.xyplot.mdguides,
        label: "show middle guides",
    });
    paper.fill = paper.checkBox({
        x: stx + 1,
        y: sty + 455,
        isChecked: commobj.xyplot.fill,
        label: "fill",
    });
    paper.jitter = paper.checkBox({
        x: stx + 50,
        y: sty + 455,
        isChecked: commobj.xyplot.jitter,
        label: "jitter points",
    });
    paper.labels = paper.checkBox({
        x: stx + 1,
        y: sty + 480,
        isChecked: commobj.xyplot.labels,
        label: "show case labels",
    });
    paper.sufnec.cover[0].click(function() {
        commobj.xyplot.sufnec = "sufficiency";
        if (xyplotdata.length > 0) {
            paper.incl.attr({"text": ("Inclusion: " + xyplotdata[3][paper.index][0])});
            paper.cov.attr({"text": ("Coverage: " + xyplotdata[3][paper.index][1])});
            paper.PRI.attr({"text": ("PRI: " + xyplotdata[3][paper.index][2])});
            paper.measures.show();
            paper.ron.hide();
        }
    });
    paper.sufnec.cover[1].click(function() {
        commobj.xyplot.sufnec = "necessity";
        if (xyplotdata.length > 0) {
            paper.incl.attr({"text": ("Inclusion: " + xyplotdata[4][paper.index][0])});
            paper.cov.attr({"text": ("Coverage: " + xyplotdata[4][paper.index][1])});
            paper.ron.attr({"text": ("Relevance: " + xyplotdata[4][paper.index][2])});
            paper.measures.show();
            paper.PRI.hide();
        }
    });
    paper.measures = paper.set();
    paper.incl = sat(paper.text(paper.sx + 2, 10, "Inclusion: "));
    paper.cov = sat(paper.text(paper.sx + 122, 10, "Coverage: "));
    paper.PRI = sat(paper.text(paper.sx + 250, 10, "PRI: "));
    paper.ron = sat(paper.text(paper.sx + 250, 10, "Relevance: "));
    paper.measures.push(paper.incl, paper.cov, paper.PRI, paper.ron);
    paper.measures.hide();
    if (commobj.xyplot.pof && xyplotdata.length > 0) {
        paper.measures.show();
        if (commobj.xyplot.sufnec == "sufficiency") {
            paper.ron.hide();
        }
    }
    paper.mdguides.cover.click(function() {
        commobj.xyplot.mdguides = paper.mdguides.isChecked;
        if (commobj.xyplot.mdguides) {
            paper.mdlines.show();
        }
        else {
            paper.mdlines.hide();
        }
    });
    paper.jitter.cover.click(function() {
        commobj.xyplot.jitter = this.isChecked;
        if (paper.xyplotdata.length > 0) {
            paper.randomjitter.x = new Array(paper.xyplotdata[0].length);
            paper.randomjitter.y = new Array(paper.xyplotdata[0].length);
            if (this.isChecked) {
                for (var i = 0; i < paper.xyplotdata[0].length; i++) {
                    paper.randomjitter.x[i] = randomBetween(-5, 5);
                    paper.randomjitter.y[i] = randomBetween(-5, 5);
                }
            }
            else {
                for (var i = 0; i < paper.xyplotdata[0].length; i++) {
                    paper.randomjitter.x[i] = 0;
                    paper.randomjitter.y[i] = 0;
                }
            }
        }
        scaleplot(paper);
        createLabels(paper);
    });
    paper.fill.cover.click(function() {
        commobj.xyplot.fill = paper.fill.isChecked;
        if (commobj.xyplot.fill) {
            paper.pointsset.attr({"fill-opacity": 1});
        }
        else {
            paper.pointsset.attr({"fill-opacity": 0});
        }
    });    
    paper.pof.cover.click(function() {
        commobj.xyplot.pof = this.isChecked;
        if (commobj.xyplot.pof && xyplotdata.length > 0) {
            paper.measures.show();
            if (commobj.xyplot.sufnec == "sufficiency") {
                paper.ron.hide();
            }
        }
        else {
            paper.measures.hide();
        }
    });
    paper.labels.cover.click(function() {
        commobj.xyplot.labels = paper.labels.isChecked;
        if (commobj.xyplot.labels) {
            paper.labelsset.show();
            paper.thsetter.show();
        }
        else {
            paper.labelsset.hide();
            paper.thsetter.hide();
        }
    });
    paper.thsetter = paper.set();
    paper.thsetter.push(sat(paper.text(stx + 20, sty + 510, "rotate")));
    paper.thsetter.push(paper.path("M" + (stx + 70) + "," + (sty + 510) + "L" + (stx + 145) + "," + (sty + 510)));
    paper.th = paper.path("M" + (stx + 70 + paper.labelRotation) + "," + (sty + 510) + "L" + (stx + 70 + paper.labelRotation - 5) + "," + (sty + 517) + "L" + (stx + 70 + paper.labelRotation + 5) + "," + (sty + 517) + "L" + (stx + 70 + paper.labelRotation) + "," + (sty + 510)).attr({"stroke-width": 1.5, fill: "#cb2626", stroke: "#cb2626"});
    paper.th.min = 0; 
    paper.th.max = 45; 
    paper.th.left = stx + 70;
    paper.th.right = stx + 145;
    paper.th.id = "xyplot";
    paper.th.drag(dragMove(paper.th), dragStart(paper.th), dragStop(paper.th));
    paper.thsetter.push(paper.th);
    paper.labelsset = paper.set();
    if (!commobj.xyplot.labels) {
        paper.labelsset.hide();
        paper.thsetter.hide();
    }
    paper.labelsArray = new Array();
    paper.x = commobj.xyplot.x;
    paper.y = commobj.xyplot.y;
    scaleplot(paper);
    createLabels(paper);
} 
} 
function draw_venn(paper) {
    var glow, txt, txtfundal;
    function hoverVenn_IN() {
        if (this.txt != "" && paper.hover) {
            glow = this.glow({
                color: "#0000ff",
                width: 2
            });
            glow.toFront();
            var BBox = this.getBBox();
            var xcoord = BBox.x;
            var ycoord = BBox.y - 20;
            if (ycoord < 0) {
                xcoord = BBox.x + 20;
                ycoord = BBox.y + 20;
            }
            txt = papers["venn"]["main"].paragraph({
                x: xcoord,
                y: ycoord,
                maxWidth: 200,
                text: this.txt.split(",").join(", "),
                textStyle: {
                    "font-family" : "Arial",
                    "font-size" : 14,
                    "text-anchor" : "start",
                    "font-weight" : "bold"
                }
            });
            var BBox2 = txt.getBBox();
            txtfundal = papers["venn"]["main"].rect(xcoord - BBox2.width/2, ycoord - 1, BBox2.width + 10, BBox2.height + 5);
            txtfundal.attr({fill: "#c9c9c9", "fill-opacity": 0.8, stroke: "none"});
            txt.toFront();
            txtfundal.translate(BBox2.width/2 - 5, -10);
            txt.show();
        }
    }
    function hoverVenn_OUT() {
        if (this.txt != "" && txt != void 0) {
            glow.remove();
            txt.remove();
            txtfundal.remove();
        }
    }
    function getCentroid(path) {
        var x = new Array(11);
        var y = new Array(11);
        var asum = 0, cxsum = 0, cysum = 0;
        var totlength = path.getTotalLength();
        for (var i = 0; i < 11; i++) {
            var location = path.getPointAtLength(i*totlength/10);
            x[i] = location.x;
            y[i] = location.y;
            if (i > 0) {
                asum += x[i - 1]*y[i] - x[i]*y[i - 1];
                cxsum += (x[i - 1] + x[i])*(x[i - 1]*y[i] - x[i]*y[i - 1]);
                cysum += (y[i - 1] + y[i])*(x[i - 1]*y[i] - x[i]*y[i - 1]);
            }
        }
        return({x: (1/(3*asum))*cxsum, y: (1/(3*asum))*cysum});
    }
    if ($("#venn").length) {
        paper.hover = true;
        if (getKeys(ttfromR).length > 0) {
            var vcolors = {
                "0": "#ffd885", 
                "1": "#96bc72", 
                "C": "#1c8ac9", 
                "?": "white"
            };
            if (paper.scale === undefined) {
                paper.scale = (Math.min($(paper.canvas).width() - 20, $(paper.canvas).height() - 70))/1000;;
            }
            if (paper.customtext === undefined) {
                paper.customtext = "";
            }
            if (paper.custom === undefined) {
                paper.custom = false;
            }
            paper.clear();
            var vennumber = ttfromR.options.conditions.length;
            var tosplit, BBox, BBox2, inIndexes, glow;
            var allshapes = paper.set();
            var boxset = paper.set();
            var borderset = paper.set();
            var labelsGroup = paper.set();
            var hoverGroup = paper.set();
            var colored = paper.set();
            var rule = paper.set();
            var customSet = paper.set();
            var customHover = paper.set();
            var customNonHover = paper.set();
            var temp = paper.rect(0, 0, 1000*paper.scale, 1000*paper.scale);
            borderset.push(temp);
            for (var i = 0; i < venn["s" + vennumber][0].length; i++) {
                var path = "M";
                var y = venn["s" + vennumber][0][i];
                for (var j = 0; j < y.length/2; j++) {
                    path += ((j == 1)?" C ":" ") + (y[2*j]*paper.scale) + "," + (y[2*j + 1]*paper.scale);
                }
                var temp = paper.path(path);
                borderset.push(temp);
            }
            for (var i = 0; i < vennumber; i++) {
                var templabel = sat(paper.text(venn[("l" + vennumber)].x[i]*paper.scale, venn[("l" + vennumber)].y[i]*paper.scale, ttfromR.options.conditions[i]), {anchor: "middle"});
                labelsGroup.push(templabel);
            }
            for (var i = 0; i < ttfromR.id.length; i++) {
                var tempath = getShape([venn["s" + vennumber][1][i]], venn["s" + vennumber][0], paper.scale);
                if (i == 0) {
                    tempath = "M 0,0 0," + 1000*paper.scale + " " + 1000*paper.scale + "," + 1000*paper.scale + " " + 1000*paper.scale + ",0 0,0 z " + tempath;
                }
                var temp = paper.path(tempath).attr({fill: vcolors[ttfromR.tt.OUT[i]], stroke: "none"});
                BBox = temp.getBBox();
                var centroid = getCentroid(temp);
                if (ttfromR.id.length == 16) {
                    if (i == 1) {
                        centroid.x += 20*paper.scale;
                        centroid.y -= 75*paper.scale;
                    }
                    else if (i == 8) {
                        centroid.x -= 20*paper.scale;
                        centroid.y -= 75*paper.scale;
                    }
                }
                var templabel;
                if (i == 0) {
                    templabel = sat(paper.text(20, 20, ttfromR.id[i]), {size: 8, anchor: "middle"});
                }
                else {
                    templabel = sat(paper.text(centroid.x, centroid.y, ttfromR.id[i]), {size: 8, anchor: "middle"});
                }
                labelsGroup.push(templabel);
                allshapes.push(temp);
                if (ttfromR.cases[i] != "") {
                    var hoverPath = temp.clone().attr({fill: "#fff", "fill-opacity": 0, stroke: "none"});
                    hoverGroup.push(hoverPath);
                    hoverPath.txt = ttfromR.tt.cases[i];
                    hoverPath.hover(hoverVenn_IN, hoverVenn_OUT, hoverPath, hoverPath);
                    colored.push(temp);
                }
            }
            allshapes.push(labelsGroup, hoverGroup, borderset);
            BBox = allshapes.getBBox();
            allshapes.transform("t10, 35");
            BBox = allshapes.getBBox();
            var colorsCols = getKeys(vcolors);
            for (var i = 0; i < colorsCols.length; i++) {
                colored.push(paper.rect(BBox.x + 50*i, BBox.y + BBox.height + 12, 11, 11)
                    .attr({fill: vcolors[colorsCols[i]]}));
                colored.push(sat(paper.text(BBox.x + 50*i + 18, BBox.y + BBox.height + 18, colorsCols[i])));
            }
            var custom = paper.checkBox({
                x: 10,
                y: 10,
                isChecked: paper.custom,
                label: "custom",
            });
            custom.cover.click(function() {
                if (custom.isChecked) {
                    paper.custom = true;
                    colored.hide();
                    hoverGroup.hide();
                    rule.show();
                    customSet.show();
                    if (glow !== undefined) {
                        glow.show();
                    }
                    customHover.toFront();
                    customNonHover.toFront();
                }
                else {
                    paper.custom = false;
                    colored.show();
                    hoverGroup.show();
                    rule.hide();
                    customSet.hide();
                    if (glow !== undefined) {
                        glow.hide();
                    }
                    hoverGroup.toFront();
                }
            });
            var ruletext = sat(paper.text(90, 16, paper.customtext), {"clip": "85, 7, 334, 20"});
            var rulerect = paper.rect(85, 6, 335, 20).attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            rule.push(ruletext, rulerect);
            if (paper.customtext != "") {
                var parsedText = parseText(paper.customtext, ttfromR.options.conditions);
                if (parsedText != "error") {
                    if (glow !== undefined) {
                        glow.remove();
                    }
                    var cols = getKeys(parsedText);
                    var inverted = new Array(cols.length);   
                    var Hovers = new Array(cols.length);
                    for (var i = 0; i < cols.length; i++) {
                        var tempinv;
                        var temp = customShape(parsedText[cols[i]],
                                               venn["s" + vennumber],
                                               paper.scale,
                                               ttfromR.id);
                        if (temp[1]) { 
                            temp[0] = paper.path("M 0,0 0," + 1000*paper.scale + " " + 1000*paper.scale + "," + 1000*paper.scale + " " + 1000*paper.scale + ",0 0,0 z " + temp[0])
                            .attr({fill: vcolors["1"], stroke: "none", "fill-opacity": 0.75});
                        }
                        else {
                            temp[0] = paper.path(temp[0]).attr({fill: vcolors["1"], stroke: "none", "fill-opacity": 0.75});
                        }
                        inverted[i] = temp[1];
                        Hovers[i] = temp[0].clone().attr({fill: "#fff", "fill-opacity": 0, stroke: "none"});
                        Hovers[i].txt = cols[i];
                        Hovers[i].hover(hoverVenn_IN, hoverVenn_OUT, Hovers[i], Hovers[i]);
                        customHover.push(Hovers[i]);
                        customSet.push(temp[0], Hovers[i]);
                    }
                    customSet.transform("t10, 35");
                    labelsGroup.toFront();
                    borderset.toFront();
                    customHover.toFront();
                    for (i = 0; i < cols.length; i++) {
                        if (!inverted[i]) {
                            Hovers[i].toFront();
                        }
                    }
                }
                else {
                    if (glow === undefined) {
                        glow = rulerect.glow({
                            color: "#ff0000",
                            width: 4
                        });
                    }
                }
                rulerect.toFront();
            }
            rulerect.click(function(e) {
                e.stopPropagation();
                var me = this;
                var BBox = this.getBBox();
                input = ruletext.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    ruletext.inlineTextEditing.stopEditing(tasta);
                    if (ruletext.attr("text") != paper.customtext) {
                        paper.customtext = ruletext.attr("text");
                        customSet.remove();
                        customHover.remove();
                        customSet = paper.set();
                        customHover = paper.set();
                        if (glow !== undefined) {
                            glow.remove();
                        }
                        if (ruletext.attr("text") != "") {
                            var tempinv;
                            var parsedText = parseText(paper.customtext, ttfromR.options.conditions);
                            if (parsedText != "error") {
                                var cols = getKeys(parsedText);
                                var inverted = new Array(cols.length);   
                                var Hovers = new Array(cols.length);
                                for (var i = 0; i < cols.length; i++) {
                                    var temp = customShape(parsedText[cols[i]],
                                               venn["s" + vennumber],
                                               paper.scale,
                                               ttfromR.id);
                                    if (temp[1]) { 
                                        temp[0] = paper.path("M 0,0 0," + 1000*paper.scale + " " + 1000*paper.scale + "," + 1000*paper.scale + " " + 1000*paper.scale + ",0 0,0 z " + temp[0])
                                        .attr({fill: vcolors["1"], stroke: "none", "fill-opacity": 0.75});
                                    }
                                    else {
                                        temp[0] = paper.path(temp[0]).attr({fill: vcolors["1"], stroke: "none", "fill-opacity": 0.75});
                                    }
                                    inverted[i] = temp[1];
                                    Hovers[i] = temp[0].clone().attr({fill: "#fff", "fill-opacity": 0}); 
                                    Hovers[i].txt = cols[i];
                                    Hovers[i].hover(hoverVenn_IN, hoverVenn_OUT, Hovers[i], Hovers[i]);
                                    customHover.push(Hovers[i]);
                                    customSet.push(temp[0], Hovers[i]);
                                }
                                customSet.transform("t10,35");
                                labelsGroup.toFront();
                                borderset.toFront();
                                customHover.toFront();
                                for (i = 0; i < cols.length; i++) {
                                    if (!inverted[i]) {
                                        Hovers[i].toFront();
                                    }
                                }
                            }
                            else {
                                glow = rulerect.glow({
                                    color: "#ff0000",
                                    width: 4
                                });
                            }
                        }
                    }
                    tasta = "enter";
                })
            });
            paper.inlineTextEditing(ruletext);
            if (paper.custom) {
                colored.hide();
                hoverGroup.hide();
                customSet.show();
                if (glow !== undefined) {
                    glow.show();
                }
            }
            else {
                colored.show();
                hoverGroup.show();
                customSet.hide();
                rule.hide();
                if (glow !== undefined) {
                    glow.hide();
                }
            }
            labelsGroup.toFront();
            borderset.toFront();
            hoverGroup.toFront();
        } 
    } 
} 
function draw_tt(paper) {
    if ($("#tt").length) {
        paper.clear();
        sat(paper.text(19, 17, "Dataset:"));
        sat(paper.text(167, 17, "Outcome:"));
        sat(paper.text(315, 17, "Conditions:"));
        var stx = 14;
        var sty = 175;
        var neg_out = paper.checkBox({
            x: stx,
            y: sty,
            isChecked: commobj.tt.neg_out,
            label: "negate outcome",
        });
        neg_out.cover.click(function() {
            commobj.tt.neg_out = neg_out.isChecked;
            console_command("tt");
        });
        var complete = paper.checkBox({
            x: stx,
            y: sty + 25,
            isChecked: commobj.tt.complete,
            label: "complete",
        });
        complete.cover.click(function() {
            commobj.tt.complete = complete.isChecked;
            console_command("tt");
        });
        var use_letters = paper.checkBox({
            x: stx,
            y: sty + 50,
            isChecked: commobj.tt.use_letters,
            label: "use letters",
        });
        use_letters.cover.click(function() {
            commobj.tt.use_letters = use_letters.isChecked;
            console_command("tt");
        });
        var show_cases = paper.checkBox({
            x: stx,
            y: sty + 75,
            isChecked: commobj.tt.show_cases,
            label: "show cases",
        });
        show_cases.cover.click(function() {
            commobj.tt.show_cases = show_cases.isChecked;
            if (commobj.tt.show_cases) {
                dcc.activate();
            }
            else {
                dcc.deactivate();
            }
            console_command("tt");
        });
        var dcc = paper.checkBox({
            x: stx,
            y: sty + 100,
            isChecked: commobj.tt.dcc,
            label: "deviant cases",
        });
        dcc.cover.click(function() {
            commobj.tt.dcc = dcc.isChecked;
            console_command("tt");
        });
        if (!commobj.tt.show_cases) {
            dcc.deactivate();
        }
       sat(paper.text(stx + 90, sty + 138, "Exclude:"), {anchor: "start"});
       paper.exclude = textbox(paper, {x: stx + 156, y: sty + 138, width: 138, height: 20, text: commobj.tt.exclude});
       paper.inlineTextEditing(paper.exclude.text);
       paper.exclude.rect.click(function(e) {
               var me = this;
               e.stopPropagation();
               var BBox = this.getBBox();
               input = paper.exclude.text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
               input.addEventListener("blur", function(e) {
                   paper.exclude.text.inlineTextEditing.stopEditing(tasta);
                   commobj.tt.exclude = paper.exclude.text.attr("text");
                   me.toFront();
                   tasta = "enter";
                   console_command("tt");
               });
       });
        sat(paper.text(stx + 160, sty + 5, "Sort by:"));
        paper.decr = sat(paper.text(stx + 240, sty + 6, "Decr."));
        paper.decr.hide();
        paper.rect(stx + 152, sty + 17.5, 78, 77)
           .attr({stroke: '#d0d0d0', 'stroke-width': 1, fill: "#ffffff", "fill-opacity": 0});
        paper.decrease = new Array(3);
        paper.rects = new Array(6);
        paper.texts = new Array(3);
        paper.positions = new Array(3);
        paper.coordsy = new Array(3);
        paper.sortsets = new Array(3);
        var keys = getKeys(commobj.tt.sort_by);
        var sortbyoptions = {"out": "outcome", "incl": "inclusion", "n": "frequency"};
        for (var i = 0; i < 3; i++) {
            paper.sortsets[i] = paper.set();
            paper.positions[i] = i;
            paper.coordsy[i] = sty + 19 + i*25;
            paper.rects[i] = paper.rect(stx + 154, paper.coordsy[i], 74, 24);
            paper.rects[i].backcolor = commobj.tt.sort_sel[keys[i]];
            paper.texts[i] = sat(paper.text(stx + 160, sty + 30 + i*25, sortbyoptions[keys[i]]));
            if (commobj.tt.sort_sel[keys[i]]) {
                paper.rects[i].attr({fill: "#79a74c", stroke: "none"});
                paper.texts[i].attr({fill: "white", "text-anchor": "start", "font-size": "14px"});
            }
            else {
                paper.rects[i].attr({fill: "#eeeeee", stroke: "none"});
                paper.texts[i].attr({fill: "black", "text-anchor": "start", "font-size": "14px"});
            }
            paper.sortsets[i].push(paper.rects[i], paper.texts[i]);
        }
        for (var i = 0; i < 3; i++) {
            paper.rects[3 + i] = paper.rect(stx + 154, paper.coordsy[i], 74, 24)
                .attr({stroke: '#ffffff', fill: "#ffffff", "fill-opacity": 0});
            paper.rects[3 + i].id = i;
            paper.rects[3 + i].name = keys[i];
            paper.rects[3 + i].top = stx + 200 - 4;
            paper.rects[3 + i].bottom = stx + 250 - 4;
            paper.sortsets[i].push(paper.rects[3 + i]);
            paper.decrease[i] = paper.checkBox({
                x: stx + 248,
                y: paper.coordsy[i] + 5,
                isChecked: commobj.tt.sort_by[keys[i]],
                label: "",
            });
            paper.decrease[i].cover.name = keys[i];
            paper.decrease[i].cover.click(function() {
                commobj.tt.sort_by[this.name] = this.isChecked;
                console_command("tt");
            });
            if (!commobj.tt.sort_sel[keys[i]]) {
                paper.decrease[i].hideIt();
            }
            paper.sortsets[i].drag(dragSortMove(paper.sortsets[i]), dragSortStart(paper.sortsets[i]), dragSortStop(paper.sortsets[i]));
        }
        if (getTrueKeys(commobj.tt.sort_sel).length == 0) {
            papers["tt"]["main"].decr.hide();
        }
        else {
            papers["tt"]["main"].decr.show();
        }
        var ctx = 396; 
        var cty = 180; 
        paper.text(ctx, cty, "Cut-off:").attr({"text-anchor": "start", "font-size": "14px"});
        paper.text(ctx - 15, cty + 25, "inclusion 1").attr({"text-anchor": "end", "font-size": "14px"});
        var ic1 = paper.text(ctx + 5, cty + 25, commobj.tt.ic1).attr({"text-anchor": "start", "font-size": "14px"});
        paper.rect(ctx, cty + 15, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(e) {
                e.stopPropagation();
                var me = this;
                var BBox = this.getBBox();
                input = ic1.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    ic1.inlineTextEditing.stopEditing(tasta);
                    if (isNumeric(ic1.attr("text"))) {
                        if (ic1.attr("text") <= 1 && ic1.attr("text") >= 0) {
                            commobj.tt.ic1 = ic1.attr("text");
                        }
                        else {
                            ic1.attr({"text": commobj.tt.ic1});
                        }
                    }
                    else {
                        ic1.attr({"text": commobj.tt.ic1});
                    }
                    me.toFront();
                    tasta = "enter";
                    console_command("tt");
                });
            });
        paper.inlineTextEditing(ic1);
        paper.text(ctx - 15, cty + 50, "inclusion 0").attr({"text-anchor": "end", "font-size": "14px"});
        var ic0 = paper.text(ctx + 5, cty + 50, commobj.tt.ic0).attr({"text-anchor": "start", "font-size": "14px"});
        paper.rect(ctx, cty + 40, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(e) {
                e.stopPropagation();
                var me = this;
                var BBox = this.getBBox();
                input = ic0.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    ic0.inlineTextEditing.stopEditing(tasta);
                    if (isNumeric(ic0.attr("text"))) {
                        if (ic0.attr("text") <= 1 && ic0.attr("text") >= 0) {
                            commobj.tt.ic0 = ic0.attr("text");
                        }
                        else {
                            ic0.attr({"text": commobj.tt.ic0});
                        }
                    }
                    else {
                        if (ic0.attr("text") == "") {
                            commobj.tt.ic0 = ic0.attr("text");
                        }
                        else {
                            ic0.attr({"text": commobj.tt.ic0});
                        }
                    }
                    me.toFront();
                    tasta = "enter";
                    console_command("tt");
                });
            });
        paper.inlineTextEditing(ic0);
        paper.text(ctx - 15, cty + 75, "PRI").attr({"text-anchor": "end", "font-size": "14px"});
        var pri_cut = paper.text(ctx + 5, cty + 75, commobj.tt.pri_cut).attr({"text-anchor": "start", "font-size": "14px"});
        paper.rect(ctx, cty + 65, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(e) {
                e.stopPropagation();
                var me = this;
                var BBox = this.getBBox();
                input = pri_cut.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    pri_cut.inlineTextEditing.stopEditing(tasta);
                    commobj.tt.pri_cut = pri_cut.attr("text");
                    me.toFront();
                    tasta = "enter";
                    console_command("tt");
                });
            });
        paper.inlineTextEditing(pri_cut);
        paper.text(ctx - 15, cty + 100, "frequency").attr({"text-anchor": "end", "font-size": "14px"});
        var frequency = paper.text(ctx + 5, cty + 100, commobj.tt.n_cut).attr({"text-anchor": "start", "font-size": "14px"});
        paper.rect(ctx, cty + 90, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(e) {
                e.stopPropagation();
                var me = this;
                var BBox = this.getBBox();
                input = frequency.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    frequency.inlineTextEditing.stopEditing(tasta);
                    commobj.tt.n_cut = frequency.attr("text");
                    me.toFront();
                    tasta = "enter";
                    console_command("tt");
                });
            });
        paper.inlineTextEditing(frequency);
        var stx = 13, sty = 0;
        var objname = paper.checkBox({
            x: stx,
            y: sty + 334 + 15,
            isChecked: commobj.tt.nameit,
            label: "Assign",
        });
        objname.cover.click(function() {
            if (this.isChecked) {
                objname.label[0].attr({"text": "Assign to:"});
                commobj.tt.objname = objname_TB.text.attr("text");
                commobj.tt.nameit = true;
                objnameset.show();
            }
            else {
                objname.label[0].attr({"text": "Assign"});
                commobj.tt.objname = "";
                commobj.tt.nameit = false;
                objnameset.hide();
            }
            console_command("tt");
        });
        var objnameset = paper.set();
        var objname_TB = textbox(paper, {x: stx + 95, y: sty + 339 + 15, width: 200, height: 20, text: commobj.tt.objname});
        paper.inlineTextEditing(objname_TB.text);
        objname_TB.rect.click(function(e) {
            var me = this;
            e.stopPropagation();
            var temp = objname_TB.text.attr("text");
            ovBox = this.getBBox();
            input = objname_TB.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function() {
                objname_TB.text.inlineTextEditing.stopEditing(tasta);
                if (objname_TB.text.attr("text") != temp) {
                    commobj.tt.objname = objname_TB.text.attr("text").replace(/[^A-Za-z0-9]/g, '');
                    if (isNumeric(commobj.tt.objname[0])) {
                        commobj.tt.objname = "x" + commobj.tt.objname;
                    }
                    objname_TB.text.attr({"text": commobj.tt.objname});
                    console_command("tt");
                }
                me.toFront();
                tasta = "enter";
            }, true);
        });
        objnameset.push(objname_TB.text, objname_TB.rect);
        if (!objname.isChecked) {
            objnameset.hide();
        }
        paper.text(ctx + 2, cty + 159 + 15, "Run").attr({"text-anchor": "start", "font-size": "14px"});
        paper.rect(ctx - 20, cty + 146 + 15, 70, 25)
        .attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
        .click(function() {
            if (info["data"][commobj["tt"].dataset].rownames != "" && !spinnerstarted) {
                console_command("tt");
                Rcommand.command = string_command;
                talkToR();
            }
        });
    }
}
function draw_minimize(paper) {
    if ($("#minimize").length) {
        paper.clear();
        var copydataset = "";
        var copyclicks = null;
        var backup = function() {
            var tempdataset = commobj["minimize"].dataset;
            commobj["minimize"].dataset = copydataset;
            copydataset = tempdataset;
            var tempclicks = null;
            if (colclicks["minimize"] !== void 0) { 
                tempclicks = copy(colclicks["minimize"]);
            }
            colclicks = copy(colclicks, exclude = "minimize");
            if (copyclicks !== null) {
                colclicks["minimize"] = copy(copyclicks);
            }
            if (tempclicks !== null) {
                copyclicks = copy(tempclicks);
            }
        }
        var dsource = paper.radio({
            x: 30,
            y: 18,
            whichChecked: 1*(commobj.minimize.source == "tt"),
            labels: ["Dataset", "TT"],
            vertspace: [0, 0],
            horspace: [0, 80],
        });
        dsource.cover[0].click(function() { 
            commobj.minimize.source = "data";
            backup();
            refresh_cols("minimize");
            paper.neg_out.activate();
            paper.use_letters.activate();
            paper.show_cases.activate();
            for (var i = 0; i < 4; i++) {
                paper.shadecutoffs[i].attr({"fill-opacity": "0"});
            }
            checkeqtt();
        });
        dsource.cover[1].click(function() { 
            commobj.minimize.source = "tt";
            backup();
            refresh_cols("minimize");
            paper.neg_out.deactivate();
            paper.use_letters.deactivate();
            paper.show_cases.deactivate();
            for (var i = 0; i < 4; i++) {
                paper.shadecutoffs[i].attr({"fill-opacity": "0.2"});
            }
            checkeqtt();
        });
        sat(paper.text(172, 17, "Outcome:"));
        sat(paper.text(320, 17, "Conditions:"));
        var stx = 14;
        var sty = 177;
        var cofftty = sty - 7;
        var txtheight = getTextHeight("0", 10);
        paper.shadecutoffs = new Array(4);
        sat(paper.text(stx + 77, cofftty + 2, "Cut-off:"));
        sat(paper.text(stx + 67, cofftty + 25, "inclusion 1"), {anchor: "end"});
        paper.ic1 = sat(paper.text(stx + 82, cofftty + 25, commobj.minimize.ic1));
        paper.shadecutoffs[0] = paper.rect(stx + 77, cofftty + 15, 50, 20, 3).attr({fill: "#000000", stroke: "none", "fill-opacity": "0.2"});
        paper.rect(stx + 77, cofftty + 15, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(event) {
                if (commobj.minimize.source == "data") {
                    event.stopPropagation();
                    var me = this;
                    var BBox = this.getBBox();
                    input = paper.ic1.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                    input.addEventListener("blur", function(e) {
                        paper.ic1.inlineTextEditing.stopEditing(tasta);
                        if (isNumeric(paper.ic1.attr("text"))) {
                            commobj.minimize.ic1 = paper.ic1.attr("text");
                            if (commobj.minimize.ic1 < commobj.minimize.ic0) {
                                commobj.minimize.ic0 = commobj.minimize.ic1;
                                paper.ic0.attr({"text": commobj.minimize.ic0});
                            }
                        }
                        else {
                            commobj.minimize.ic1 = "1";
                            commobj.minimize.ic0 = "";
                            paper.ic1.attr({"text": "1"});
                            paper.ic0.attr({"text": ""});
                        }
                        me.toFront();
                        tasta = "enter";
                        console_command("minimize");
                    });
                }
            });
        paper.inlineTextEditing(paper.ic1);
        sat(paper.text(stx + 67, cofftty + 50, "inclusion 0"), {anchor: "end"});
        paper.ic0 = sat(paper.text(stx + 82, cofftty + 50, commobj.minimize.ic0));
        paper.shadecutoffs[1] = paper.rect(stx + 77, cofftty + 40, 50, 20, 3).attr({fill: "#000000", stroke: "none", "fill-opacity": "0.2"});
        paper.rect(stx + 77, cofftty + 40, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(event) {
                if (commobj.minimize.source == "data") {
                    event.stopPropagation();
                    var me = this;
                    var BBox = this.getBBox();
                    input = paper.ic0.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                    input.addEventListener("blur", function(e) {
                        paper.ic0.inlineTextEditing.stopEditing(tasta);
                        commobj.minimize.ic0 = paper.ic0.attr("text");
                        if (isNumeric(paper.ic0.attr("text"))) {
                            commobj.minimize.ic0 = paper.ic0.attr("text");
                            if (commobj.minimize.ic1 < commobj.minimize.ic0) {
                                commobj.minimize.ic0 = commobj.minimize.ic1;
                                paper.ic0.attr({"text": commobj.minimize.ic0});
                            }
                        }
                        else {
                            commobj.minimize.ic0 = commobj.minimize.ic1;
                            paper.ic0.attr({"text": commobj.minimize.ic0});
                        }
                        me.toFront();
                        tasta = "enter";
                        console_command("minimize");
                    });
                }
            });
        paper.inlineTextEditing(paper.ic0);
        sat(paper.text(stx + 67, cofftty + 75, "PRI"), {anchor: "end"});
        paper.pri_cut = sat(paper.text(stx + 82, cofftty + 75, commobj.minimize.pri_cut));
        paper.shadecutoffs[2] = paper.rect(stx + 77, cofftty + 65, 50, 20, 3).attr({fill: "#000000", stroke: "none", "fill-opacity": "0.2"});
        paper.rect(stx + 77, cofftty + 65, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(event) {
                if (commobj.minimize.source == "data") {
                    event.stopPropagation();
                    var me = this;
                    var BBox = this.getBBox();
                    input = paper.pri_cut.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                    input.addEventListener("blur", function(e) {
                        paper.pri_cut.inlineTextEditing.stopEditing(tasta);
                        commobj.minimize.pri_cut = paper.pri_cut.attr("text");
                        me.toFront();
                        tasta = "enter";
                        console_command("minimize");
                    });
                }
            });
        paper.inlineTextEditing(paper.pri_cut);
        sat(paper.text(stx + 67, cofftty + 100, "frequency"), {anchor: "end"});
        paper.n_cut = sat(paper.text(stx + 82, cofftty + 100, commobj.minimize.n_cut));
        paper.shadecutoffs[3] = paper.rect(stx + 77, cofftty + 90, 50, 20, 3).attr({fill: "#000000", stroke: "none", "fill-opacity": "0.2"});
        paper.rect(stx + 77, cofftty + 90, 50, 20, 3)
            .attr({fill: "#ffffff", stroke: "#a0a0a0", "fill-opacity": "0"})
            .click(function(event) {
                if (commobj.minimize.source == "data") {
                    event.stopPropagation();
                    var me = this;
                    var BBox = this.getBBox();
                    input = paper.n_cut.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                    input.addEventListener("blur", function(e) {
                        paper.n_cut.inlineTextEditing.stopEditing(tasta);
                        commobj.minimize.n_cut = paper.n_cut.attr("text");
                        me.toFront();
                        tasta = "enter";
                        console_command("minimize");
                    });
                }
            });
        paper.inlineTextEditing(paper.n_cut);
        paper.neg_out = paper.checkBox({
            x: stx,
            y: sty + 124,
            isChecked: commobj.minimize.neg_out,
            label: "negate outcome",
        });
        paper.neg_out.cover.click(function() {
            if (this.active) {
                commobj.minimize.neg_out = paper.neg_out.isChecked;
                console_command("minimize");
            }
        });
        paper.show_cases = paper.checkBox({
            x: stx,
            y: sty + 149,
            isChecked: commobj.minimize.show_cases,
            label: "show cases",
        });
        paper.show_cases.cover.click(function() {
            if (this.active) {
                commobj.minimize.show_cases = paper.show_cases.isChecked;
                console_command("minimize");
            }
        });
        paper.use_letters = paper.checkBox({
            x: stx,
            y: sty + 174,
            isChecked: commobj.minimize.use_letters,
            label: "use letters",
        });
        paper.use_letters.cover.click(function() {
            if (this.active) {
                commobj.minimize.use_letters = paper.use_letters.isChecked;
                console_command("minimize");
            }
        });
        paper.neg_out.deactivate();
        paper.use_letters.deactivate();
        paper.show_cases.deactivate();
        stx = 115;
        sty = 195;
        paper.dir_exp_lbl = sat(paper.text(stx + 145, sty - 23, "")); 
        paper.dir_exp = paper.set();
        var dir_exp_txt = sat(paper.text(stx + 125, sty, commobj.minimize.dir_exp));
        paper.dir_exp.push(dir_exp_txt);
        paper.dir_exp.push(sat(paper.rect(stx + 120, sty - 10, 218, 20, 3))
            .click(function(e) {
                var me = this;
                e.stopPropagation();
                var BBox = this.getBBox();
                input = dir_exp_txt.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
                input.addEventListener("blur", function(e) {
                    dir_exp_txt.inlineTextEditing.stopEditing(tasta);
                    commobj.minimize.dir_exp = dir_exp_txt.attr("text");
                    me.toFront();
                    tasta = "enter";
                    console_command("minimize");
                });
            }));
        paper.inlineTextEditing(dir_exp_txt);
        sat(paper.text(stx + 55, sty - 23, "Include"));
        paper.rect(stx + 52, sty - 11, 50, 23)
                   .attr({stroke: '#d0d0d0', 'stroke-width': 1, fill: "#ffffff", "fill-opacity": 0});
        var rects = new Array(4);
        var texts = new Array(2);
        var selected = false;
        selected = commobj.minimize.include[0] == "?";
        rects[0] = paper.rect(stx + 54, sty - 9, 22, 19)
                   .attr({fill: selected?"#79a74c":"#eeeeee", stroke: 'none'});
        texts[0] = sat(paper.text(stx + 61, sty + 1, "?"));
        rects[0].backcolor = selected;
        if (selected) {
            texts[0].attr({"fill": "white"});
            paper.dir_exp_lbl.attr({"text": "Directional expectations:"});
            paper.dir_exp.show();
        }
        else {
            paper.dir_exp_lbl.attr({"text": ""});
            paper.dir_exp.hide();
        }
        selected = commobj.minimize.include[1] == "C";
        rects[1] = paper.rect(stx + 78, sty - 9, 22, 19)
                   .attr({fill: selected?"#79a74c":"#eeeeee", stroke: 'none'});
        texts[1] = sat(paper.text(stx + 84, sty + 1, "C"));
        rects[1].backcolor = selected;
        if (selected) {
            texts[1].attr({"fill": "white"});
        }
        rects[2] = paper.rect(stx + 53, sty - 8, 23, 21)
                   .attr({stroke: 'none', fill: "#ffffff", "fill-opacity": 0, cursor: "pointer"});
        rects[2].click(function() {
                commobj.minimize.pi_cons = "0"; 
                commobj.minimize.sol_cons = "0"; 
                paper.sol_cons.text.attr({"text": "0"}); 
                paper.pi_cons.text.attr({"text": "0"}); 
                if (commobj.minimize.include[0] == "?") {
                    rects[0].attr({fill: "#eeeeee", stroke: "none"});
                    texts[0].attr({"text-anchor": "start", "font-size": "14px", fill: "black"});
                    commobj.minimize.include[0] = "";
                    paper.dir_exp_lbl.attr({"text": ""});
                    paper.dir_exp.hide();
                }
                else {
                    rects[0].attr({fill: "#79a74c", stroke: "none"});
                    texts[0].attr({"text-anchor": "start", "font-size": "14px", fill: "white"});
                    commobj.minimize.include[0] = "?";
                    paper.dir_exp_lbl.attr({"text": "Directional expectations:"});
                    paper.dir_exp.show();
                }
                console_command("minimize");
        });
        rects[3] = paper.rect(stx + 78, sty - 8, 23, 21)
                   .attr({stroke: 'none', fill: "#ffffff", "fill-opacity": 0, cursor: "pointer"});
        rects[3].click(function() {
                if (commobj.minimize.include[1] == "C") {
                    rects[1].attr({fill: "#eeeeee", stroke: "none"});
                    texts[1].attr({"text-anchor": "start", "font-size": "14px", fill: "black"});
                    commobj.minimize.include[1] = "";
                }
                else {
                    rects[1].attr({fill: "#79a74c", stroke: "none"});
                    texts[1].attr({"text-anchor": "start", "font-size": "14px", fill: "white"});
                    commobj.minimize.include[1] = "C";
                }
                console_command("minimize");
        });
        stx = 174;
        sty = 228;
        var details = paper.checkBox({
            x: stx,
            y: sty,
            isChecked: commobj.minimize.details,
            label: "show details",
        });
        details.cover.click(function() {
            commobj.minimize.details = details.isChecked;
            console_command("minimize");
        });
        var row_dom = paper.checkBox({
            x: stx,
            y: sty + 25,
            isChecked: commobj.minimize.row_dom,
            label: "row dominance",
        });
        row_dom.cover.click(function() {
            commobj.minimize.row_dom = row_dom.isChecked;
            console_command("minimize");
        });
        var all_sol = paper.checkBox({
            x: stx + 130,
            y: sty,
            isChecked: commobj.minimize.all_sol,
            label: "maximal solutions",
        });
        all_sol.cover.click(function() {
            commobj.minimize.all_sol = all_sol.isChecked;
            console_command("minimize");
        });
        var first_min = paper.checkBox({
            x: stx + 130,
            y: sty + 25,
            isChecked: commobj.minimize.first_min,
            label: "minimal PI no.",
        });
        first_min.cover.click(function() {
            commobj.minimize.first_min = first_min.isChecked;
            console_command("minimize");
        });
        if (!commobj.minimize.details && commobj.minimize.source == "data") {
            paper.show_cases.deactivate();
        }
        stx = 240;
        sty = 189;
        sat(paper.text(stx + 60, sty + 105 + 13, "Search solutions by:"), {anchor: "end"});
        sat(paper.text(stx, sty + 130 + 13, "coverage"), {anchor: "end"});
        sat(paper.text(stx, sty + 155 + 13, "consistency"), {anchor: "end"});
        sat(paper.text(stx, sty + 180 + 13, "PI consistency"), {anchor: "end"});
        paper.sol_cov = textbox(paper, {x: stx + 15, y: sty + 130 + 13, width: 50, height: 20, text: commobj.minimize.sol_cov});
        paper.inlineTextEditing(paper.sol_cov.text);
        paper.sol_cov.rect.click(function(e) {
            var me = this;
            e.stopPropagation();
            var BBox = this.getBBox();
            input = paper.sol_cov.text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                paper.sol_cov.text.inlineTextEditing.stopEditing(tasta);
                commobj.minimize.sol_cov = paper.sol_cov.text.attr("text");
                me.toFront();
                tasta = "enter";
                console_command("minimize");
            });
        });
        paper.sol_cons = textbox(paper, {x: stx + 15, y: sty + 155 + 13, width: 50, height: 20, text: commobj.minimize.sol_cons});
        paper.inlineTextEditing(paper.sol_cons.text);
        paper.sol_cons.rect.click(function(e) {
            var me = this;
            e.stopPropagation();
            var BBox = this.getBBox();
            input = paper.sol_cons.text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                paper.sol_cons.text.inlineTextEditing.stopEditing(tasta);
                commobj.minimize.sol_cons = paper.sol_cons.text.attr("text");
                if (commobj.minimize.sol_cons != "1") {
                    rects[0].attr({fill: "#79a74c", stroke: "none"});
                    texts[0].attr({"text-anchor": "start", "font-size": "14px", fill: "white"});
                    commobj.minimize.include[0] = "?";
                }
                me.toFront();
                tasta = "enter";
                console_command("minimize");
            });
        });
        paper.pi_cons = textbox(paper, {x: stx + 15, y: sty + 180 + 13, width: 50, height: 20, text: commobj.minimize.pi_cons});
        paper.inlineTextEditing(paper.pi_cons.text);
        paper.pi_cons.rect.click(function(e) {
            var me = this;
            e.stopPropagation();
            var BBox = this.getBBox();
            input = paper.pi_cons.text.inlineTextEditing.startEditing(BBox.x + 1, BBox.y + 21 - 1*(navigator.browserType == "Firefox"), BBox.width - 2, BBox.height - 2);
            input.addEventListener("blur", function(e) {
                paper.pi_cons.text.inlineTextEditing.stopEditing(tasta);
                commobj.minimize.pi_cons = paper.pi_cons.text.attr("text");
                if (commobj.minimize.pi_cons != "1") {
                    rects[0].attr({fill: "#79a74c", stroke: "none"});
                    texts[0].attr({"text-anchor": "start", "font-size": "14px", fill: "white"});
                    commobj.minimize.include[0] = "?";
                }
                me.toFront();
                tasta = "enter";
                console_command("minimize");
            });
        });
        stx = 380;
        sty = 310;
        sat(paper.text(stx + 4, sty + 13, "Depth:"));
        paper.soldepth = paper.counter({"x": stx + 25, "y": sty + 38, "startval": 0, "maxval": 9, "width": 18, "textheight": txtheight});
        paper.soldepth.label({"label": "solution", "x": -32, "y": 0, "anchor": "end"}); 
        paper.pidepth = paper.counter({"x": stx + 25, "y": sty + 63, "startval": 0, "maxval": 9, "width": 18, "textheight": txtheight});
        paper.pidepth.label({"label": "PI", "x": -32, "y": 0, "anchor": "end"});
        paper.soldepth.update = function() {
            commobj.minimize.sol_depth = paper.soldepth.value;
            console_command("minimize");
        }
        paper.soldepth.down.click(function() {
            paper.soldepth.update();
        });
        paper.soldepth.up.click(function() {
            paper.soldepth.update();
        });
        paper.pidepth.update = function() {
            commobj.minimize.pi_depth = paper.pidepth.value;
            console_command("minimize");
        }
        paper.pidepth.down.click(function() {
            paper.pidepth.update();
        });
        paper.pidepth.up.click(function() {
            paper.pidepth.update();
        });
        stx = 14;
        sty = 180 + 25;
        var objname = paper.checkBox({
            x: stx,
            y: sty + 215,
            isChecked: commobj.minimize.nameit,
            label: "Assign",
        });
        objname.cover.click(function() {
            if (this.isChecked) {
                objname.label[0].attr({"text": "Assign to:"});
                commobj.minimize.objname = objname_TB.text.attr("text");
                commobj.minimize.nameit = true;
                objnameset.show();
            }
            else {
                objname.label[0].attr({"text": "Assign"});
                commobj.minimize.objname = "";
                commobj.minimize.nameit = false;
                objnameset.hide();
            }
            console_command("minimize");
        });
        var objnameset = paper.set();
        var objname_TB = textbox(paper, {x: stx + 91, y: sty + 221, width: 200, height: 20, text: commobj.minimize.objname});
        paper.inlineTextEditing(objname_TB.text);
        objname_TB.rect.click(function(e) {
            var me = this;
            e.stopPropagation();
            var temp = objname_TB.text.attr("text");
            ovBox = this.getBBox();
            input = objname_TB.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function() {
                objname_TB.text.inlineTextEditing.stopEditing(tasta);
                if (objname_TB.text.attr("text") != temp) {
                    commobj.minimize.objname = objname_TB.text.attr("text").replace(/[^A-Za-z0-9]/g, '');
                    if (isNumeric(commobj.minimize.objname[0])) {
                        commobj.minimize.objname = "x" + commobj.minimize.objname;
                    }
                    objname_TB.text.attr({"text": commobj.minimize.objname});
                    if (commobj.minimize.dataset != "") {
                        console_command("minimize");
                    }
                }
                me.toFront();
                tasta = "enter";
            }, true);
        });
        objnameset.push(objname_TB.text, objname_TB.rect);
        if (!objname.isChecked) {
            objnameset.hide();
        }
        sat(paper.text(stx + 383, sty + 220, "Run"));
        paper.rect(stx + 360, sty + 207, 70, 25)
        .attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
        .click(function() {
            if (commobj.minimize.source == "tt") {
                if (info["tt"] !== null) {
                    if (info["tt"][commobj.minimize.dataset] !== void 0 && !spinnerstarted) {
                        console_command("minimize");
                        Rcommand.command = string_command;
                        talkToR();
                    }
                }
            }
            else {
                if (info["data"][commobj.minimize.dataset].rownames != "" && !spinnerstarted) {
                    console_command("minimize");
                    commobj.minimize.outcome = getTrueKeys(colclicks.minimize.outcome);
                    commobj.minimize.conditions = getTrueKeys(colclicks.minimize.conditions);
                    Rcommand.command = string_command;
                    talkToR();
                }
            }
        });
    }
}
function draw_findRows(paper) {
    if ($("#findRows").length) {
        paper.clear();
        sat(paper.text(19, 17, "Truth table:"));
        var type = new Array(4);
        type[0] = paper.checkBox({
            x: 170,
            y: 30,
            isChecked: commobj.findRows.type[0],
            label: "Select All",
        });
        type[1] = paper.checkBox({
            x: 170,
            y: 55,
            isChecked: commobj.findRows.type[1],
            label: "Expression subsets:",
        });
        type[2] = paper.checkBox({
            x: 170,
            y: 110,
            isChecked: commobj.findRows.type[2],
            label: "Contradictory simplifying assumptions",
        });
        type[3] = paper.checkBox({
            x: 170,
            y: 135,
            isChecked: commobj.findRows.type[3],
            label: "Simultaneous subset relations",
        });
        var remainders = paper.checkBox({
            x: 430,
            y: 55,
            isChecked: commobj.findRows.remainders,
            label: "remainders only",
            pos: 1,
        });
        var shade = paper.rect(195 - 5, 88 - 10, 255, 20, 3).attr({fill: "#000000", stroke: "none", "fill-opacity": "0.2"});
        var type1_TB = textbox(paper, {x: 195, y: 88, width: 255, height: 20, text: commobj.findRows.expression});
        paper.inlineTextEditing(type1_TB.text);
        type1_TB.rect.click(function(e) {
            if (commobj.findRows.type[1]) { 
                var me = this;
                e.stopPropagation();
                var temp = type1_TB.text.attr("text");
                ovBox = this.getBBox();
                input = type1_TB.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
                input.addEventListener("blur", function() {
                    type1_TB.text.inlineTextEditing.stopEditing(tasta);
                    if (type1_TB.text.attr("text") != temp) {
                        commobj.findRows.expression = type1_TB.text.attr("text");
                        type1_TB.text.attr({"text": commobj.findRows.expression});
                        console_command("findRows");
                    }
                    me.toFront();
                    tasta = "enter";
                }, true);
            }
        });
        if (type[1].isChecked) {
            shade.attr({"fill-opacity": "0"});
        }
        type[0].cover.click(function() {
            var isChecked = type[0].isChecked;
            shade.attr({"fill-opacity": 1*(isChecked ? 0 : 0.2)});
            for (var i = 0; i < 4; i++) {
                commobj.findRows.type[i] = isChecked;
                if (isChecked) {
                    type[i].check();
                    remainders.activate();
                }
                else {
                    type[i].uncheck();
                    remainders.deactivate();
                }
            }
            console_command("findRows");
        });
        var allchecked = function() {
            var check = true;
            for (var i = 1; i < 4; i++) {
                check = check && type[i].isChecked;
            }
            return(check);
        }
        type[1].cover.click(function() {
            commobj.findRows.type[1] = type[1].isChecked;
            shade.attr({"fill-opacity": 1*(type[1].isChecked ? 0 : 0.2)});
            if (allchecked()) {
                type[0].check();
                commobj.findRows.type[0] = true;
                shade.attr({"fill-opacity": 0});
            }
            else {
                type[0].uncheck();
                commobj.findRows.type[0] = false;
            }
            if (type[1].isChecked) {
                remainders.activate();
            }
            else {
                remainders.deactivate();
            }
            console_command("findRows");
        });
        type[2].cover.click(function() {
            commobj.findRows.type[2] = type[2].isChecked;
            if (allchecked()) {
                type[0].check();
                commobj.findRows.type[0] = true;
            }
            else {
                type[0].uncheck();
                commobj.findRows.type[0] = false;
            }
            console_command("findRows");
        });
        type[3].cover.click(function() {
            commobj.findRows.type[3] = type[3].isChecked;
            if (allchecked()) {
                type[0].check();
                commobj.findRows.type[0] = true;
            }
            else {
                type[0].uncheck();
                commobj.findRows.type[0] = false;
            }
            console_command("findRows");
        });
        remainders.cover.click(function() {
            commobj.findRows.remainders = remainders.isChecked;
        });
        var objname = paper.checkBox({
            x: 14,
            y: 170,
            isChecked: commobj.findRows.nameit,
            label: "Assign",
        });
        objname.cover.click(function() {
            if (this.isChecked) {
                objname.label[0].attr({"text": "Assign to:"});
                commobj.findRows.objname = objname_TB.text.attr("text");
                commobj.findRows.nameit = true;
                objnameset.show();
            }
            else {
                objname.label[0].attr({"text": "Assign"});
                commobj.findRows.objname = "";
                commobj.findRows.nameit = false;
                objnameset.hide();
            }
            console_command("findRows");
        });
        var objnameset = paper.set();
        var objname_TB = textbox(paper, {x: 14 + 95, y: 170 + 6, width: 220, height: 20, text: commobj.findRows.objname});
        paper.inlineTextEditing(objname_TB.text);
        objname_TB.rect.click(function(e) {
            var me = this;
            e.stopPropagation();
            var temp = objname_TB.text.attr("text");
            ovBox = this.getBBox();
            input = objname_TB.text.inlineTextEditing.startEditing(ovBox.x + 1, ovBox.y + 21 - 1*(navigator.browserType == "Firefox"), ovBox.width - 2, ovBox.height - 2);
            input.addEventListener("blur", function() {
                objname_TB.text.inlineTextEditing.stopEditing(tasta);
                if (objname_TB.text.attr("text") != temp) {
                    commobj.findRows.objname = objname_TB.text.attr("text").replace(/[^A-Za-z0-9]/g, '');
                    if (isNumeric(commobj.findRows.objname[0])) {
                        commobj.findRows.objname = "x" + commobj.findRows.objname;
                    }
                    objname_TB.text.attr({"text": commobj.findRows.objname});
                    console_command("findRows");
                }
                me.toFront();
                tasta = "enter";
            }, true);
        });
        objnameset.push(objname_TB.text, objname_TB.rect);
        if (!objname.isChecked) {
            objnameset.hide();
        }
        sat(paper.text(377 + 23, 163 + 13, "Run"));
        paper.rect(377, 163, 70, 25)
        .attr({fill: "white", "fill-opacity": 0, 'stroke-width': 1.25})
        .click(function() {
            if (info["tt"][commobj.findRows.obj] !== void 0 && !spinnerstarted) {
                Rcommand.command = string_command;
                talkToR();
            }
        });
        console_command("findRows");
    }
}
function checkIfDataLoadedInR() {
    updatecounter += 1;
    if (updatecounter < 101) { 
        if (tempdatainfo.nrows > 0) {
            refresh_cols("import");
            if (!commobj["read_table"].customname) {
                papers["import"]["main"].objnametext.attr({"text": dirfile.filename});
                commobj["read_table"].objname = dirfile.filename;
            }
            else {
                commobj["read_table"].objname = papers["import"]["main"].objnametext.attr("text");
            }
            console_command("import");
            updatecounter = 0; 
        }
        else {
            setTimeout(checkIfDataLoadedInR, 50);
        }
    }
    else {
        if (rloadcycles == 0) {
            $("#result_main").append("<br><br><span style='color:red'>Warning: R takes rather long to load the data...</span><br>");
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1000);
            rloadcycles += 1;
            setTimeout(checkIfDataLoadedInR, 50);
        }
        else if (rloadcycles == 1) {
            $("#result_main").append("<br><span style='color:red'>Warning: R still takes long to load the data...</span><br>");
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1000);
            rloadcycles += 1;
            setTimeout(checkIfDataLoadedInR, 50);
        }
        else if (rloadcycles == 2) {
            $("#result_main").append("<br><span style='color:red'>Warning: is this such a big dataset...?</span><br>");
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1000);
            rloadcycles += 1;
            setTimeout(checkIfDataLoadedInR, 50);
        }
        else if (rloadcycles == 3) {
            $("#result_main").append("<br><span style='color:red'>Warning: we're getting bored here...</span><br>");
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1000);
            rloadcycles += 1;
            setTimeout(checkIfDataLoadedInR, 50);
        }
        else if (rloadcycles == 4) {
            $("#result_main").append("<br><span style='color:red'>Error: OK I give up. It's just too long.</span><br>");
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1000);
            rloadcycles = 0;
        }
        updatecounter = 0; 
    }
}
function consoleIfPathChanges() {
    updatecounter2 += 1;
    if (updatecounter2 < 21) { 
        if (dirfile.filepath[0][0] != pathcopy[0][0]) {
            console_command(current_command);
            updatecounter2 = 0; 
        }
        else {
            setTimeout(consoleIfPathChanges, 50);
        }
    }
    else {
        updatecounter2 = 0; 
    }
}
function doWhenDataPointsAreReturned() {
    updatecounter += 1;
    if (updatecounter < 101) { 
        if (responseR) {
            if (poinths.message == "notnumeric") {
                commobj.calibrate.thresholds = new Array(6);
                for (var i = 0; i < 6; i++) {
                    ths[i].attr({"text": ""});
                    commobj.calibrate.thresholds[i] = "";
                }
            }
            else if (poinths.message == "OK") {
                if (poinths.thvals.length > 0) {
                    for (var i = 0; i < 6; i++) {
                        ths[i].attr({"text": ""});
                        commobj.calibrate.thresholds[i] = "";
                        commobj.calibrate.thscopycrp[i] = "";
                    }
                    for (var i = 0; i < poinths.thvals.length; i++) {
                        ths[i].attr({"text": poinths.thvals[i]});
                        commobj.calibrate.thresholds[i] = poinths.thvals[i];
                        commobj.calibrate.thscopycrp[i] = poinths.thvals[i];
                    }
                }
            }
            if (focus == "calibrate") {
                console_command("calibrate");
            }
            drawPointsAndThresholds();
        }
        else {
            setTimeout(doWhenDataPointsAreReturned, 50);
        }
    }
    else {
        updatecounter = 0; 
    }
}
function drawPointsAndThresholds() {
    if (colclicks.calibrate !== void 0) {
        thsetter_content.remove();
        thsetter_content = papers["calibrate"]["main"].set().attr({stroke: "#a0a0a0"});
        if (poinths.vals.length > 0) {
            var pointsmin = min(poinths.vals);
            var pointsmax = max(poinths.vals);
            var mint = (any(isNumeric(commobj.calibrate.thresholds), "==", true) ? min(commobj.calibrate.thresholds) : pointsmin);
            var maxt = (any(isNumeric(commobj.calibrate.thresholds), "==", true) ? max(commobj.calibrate.thresholds) : pointsmax);
            var minv = min([pointsmin, poinths.prettyx[0],                          ((papers["calibrate"]["main"].crfuz == 0) ? pointsmin : min([mint, pointsmin]))]);
            var maxv = max([pointsmax, poinths.prettyx[poinths.prettyx.length - 1], ((papers["calibrate"]["main"].crfuz == 0) ? pointsmax : max([maxt, pointsmax]))]);
            var lm = 165;
            var rm = $("#calibrate").width() - 20;
            var thy = 304;
            yvals = rep(280, poinths.vals.length);
            var xs = new Array(poinths.prettyx.length);
            var patharray = new Array(2*(xs.length + 1));
            for (var i = 0; i < poinths.prettyx.length; i++) {
                xs[i] = (rm - lm)*(poinths.prettyx[i] - minv)/(maxv - minv) + lm;
                patharray[2*i] = ["M", xs[i], thy + 4.5];
                patharray[2*i + 1] = ["L", xs[i], thy - 0.5];
            }
            patharray[patharray.length - 2] = ["M", xs[0] - 0.5, thy - 0.5];
            patharray[patharray.length - 1] = ["L", xs[xs.length - 1] + 0.5, thy - 0.5];
            for (var i = 0; i < xs.length; i++) {
                thsetter_content.push(sat(papers["calibrate"]["main"].text(xs[i], thy + 15, poinths.prettyx[i]), {"size": 10, "anchor": "middle"}))
            }
            thsetter_content.push(papers["calibrate"]["main"].path(patharray)); 
            if (papers["calibrate"]["main"].crfuz == 0) {
                if (commobj.calibrate.findth && poinths.thvals !== void 0) {
                    commobj.calibrate.thresholds = ["", "", "", "", "", ""];
                    for (var i = 0; i < poinths.thvals.length; i++) {
                        commobj.calibrate.thresholds[i] = poinths.thvals[i];
                        ths[i].attr({"text": poinths.thvals[i]});
                    }
                }
                if (commobj.calibrate.jitter) {
                    for (var i = 0; i < poinths.vals.length; i++) {
                        yvals[i] = randomBetween(250, 280);
                    }
                }
            }
            else {
                if (poinths.vals.length == poinths.fuzzyvals.length) {
                    vertax = [0, 0.5, 1];
                    vertaxl = ["0", ".5", "1"];
                    patharray = new Array(2*(vertax.length + 1));
                    for (var i = 0; i < vertax.length; i++) {
                        thsetter_content.push(sat(papers["calibrate"]["main"].text(lm - 14, thy - 5 - 110*vertax[i], vertaxl[i]), {"size": 10, "anchor": "end"}));
                        patharray[2*i] = ["M", lm - 10, thy - 5 - 110*vertax[i]];
                        patharray[2*i + 1] = ["L", lm - 5, thy - 5 - 110*vertax[i]]; 
                    }
                    patharray[patharray.length - 2] = ["M", lm - 5, thy - 4.5 - 110*0];
                    patharray[patharray.length - 1] = ["L", lm - 5, thy - 5.5 - 110*1];
                    thsetter_content.push(papers["calibrate"]["main"].path(patharray)); 
                    for (var i = 0; i < poinths.vals.length; i++) {
                        yvals[i] = thy - 5 - 110*poinths.fuzzyvals[i];
                    }
                }
            }
            var position, th;
            for (var i = 0; i < 6; i++) {
                if (commobj.calibrate.thresholds[i] != "") {
                    if (commobj.calibrate.thresholds[i] != void 0) {
                        if (commobj.calibrate.thresholds[i] < minv) {
                            commobj.calibrate.thresholds[i] = minv;
                            ths[i].attr({"text": minv});
                        }
                        if (commobj.calibrate.thresholds[i] > maxv) {
                            commobj.calibrate.thresholds[i] = maxv;
                            ths[i].attr({"text": maxv});
                        }
                        position = (rm - lm)*(commobj.calibrate.thresholds[i] - minv)/(maxv - minv) + lm;
                        handles[i] = papers["calibrate"]["main"].path([
                            ["M", position, thy],
                            ["L", position - 5, thy + 7],
                            ["L", position + 5, thy + 7],
                            ["L", position, thy],
                            ["L", position, thy - 120]
                        ]).attr({"stroke-width": 1, fill: "#cb2626", stroke: "#cb2626"});
                        handles[i].min = ((papers["calibrate"]["main"].crfuz == 0) ? pointsmin : min([pointsmin, mint]));
                        handles[i].max = ((papers["calibrate"]["main"].crfuz == 0) ? pointsmax : max([pointsmax, maxt]));
                        handles[i].name = i;
                        handles[i].left = (rm - lm)*(handles[i].min - minv)/(maxv - minv) + lm;
                        handles[i].right = (rm - lm)*(handles[i].max - minv)/(maxv - minv) + lm;
                        handles[i].id = "thsetter";
                        handles[i].pos = position;
                        handles[i].drag(dragMove(handles[i]), dragStart(handles[i]), dragStop(handles[i]));
                        thsetter_content.push(handles[i]);
                    }
                }
            }
            for (var i = 0; i < poinths.vals.length; i++) {
                var point = papers["calibrate"]["main"].circle(
                    (rm - lm)*(poinths.vals[i] - minv)/(maxv - minv) + lm,
                    yvals[i], 3);
                point.attr({fill: "#ffffff", "fill-opacity": 0.0});
                point.txt = info["data"][commobj["calibrate"].dataset].rownames[i];
                point.hover(hoverIn, hoverOut, point, point);
                thsetter_content.push(point);
            }
            var txt, txtfundal;
            function hoverIn() {
                var BBox = this.getBBox();
                var xcoord = BBox.x;
                var ycoord = BBox.y - 20;
                txt = sat(papers["calibrate"]["main"].text(xcoord, ycoord, this.txt), {"anchor": "middle"});
                txt.attr({"font-weight": "bold", "fill-opacity": 0.7});
                var BBox2 = txt.getBBox();
                txtfundal = papers["calibrate"]["main"].rect(xcoord - BBox2.width/2, ycoord - 1, BBox2.width + 10, 16);
                txtfundal.attr({fill: "#c9c9c9", "fill-opacity": 0.6, stroke: "none"});
                txt.toFront();
                txt.translate(5, 7);
                txt.attr({"font-weight": "bold"});
                txt.show();
            }
            function hoverOut() {
                txt.remove();
                txtfundal.remove();
            }
        }
    }
}
function doWhenXYplotPointsAreReturned() {
    updatecounter += 1;
    if (updatecounter < 21) {
        if (lastvals.toString() != xyplotdata.toString()) {
            draw_xyplot(papers["xyplot"]["main"]);
            updatecounter = 0;
        }
        else {
            setTimeout(doWhenXYplotPointsAreReturned, 50);
        }
    }
    else {
        $("#result").append("<br><br><span style='color:red'>Error in doWhenXYplotPointsAreReturned:<br> R takes too long to respond.</span><br>");
        $("#result").animate({
            scrollTop: $("#result")[0].scrollHeight
        }, 1000);
        updatecounter = 0; 
    }
}
function pingit(where) {
    if (where == "loop") { 
        pingobj += 1;
        Shiny.onInputChange("pingobj", pingobj);
        setTimeout(pingit.bind(null, "loop"), 1000*60*4);
    }
    else {
        if (!pingstarted) {
            pingstarted = true;
            pingobj += 1;
            Shiny.onInputChange("pingobj", pingobj);
            setTimeout(pingit.bind(null, "loop"), 1000*60*4);
        }
    }
}
function printRcommand() {
    updatecounter += 1;
    if (responseR) {
        var viewdata = "";
        txtcommand = ""
        updatecounter = 0;
        for (var i = 0, length = outres.evaluate.length; i < length; i++) {
            var temp = outres.evaluate[i].command.split("\n");
            var first = ">";
            var color = "#932192";
            if (tempcommand != "") {
                temp.splice(0, tempcommand.split("\n").length - 1);
                first = "+";
                color = "blue";
            }
            for (var t = 0, tlen = temp.length; t < tlen; t++) {
                temp[t] = strwrap(temp[t], cpl.result, "  ");
            }
            $("#result_main").append("<span style='color:" + color + "'>" + first + " </span>");
            $("#result_main").append("<span style='color:blue'>" + 
                temp.join("<br>+ ").split(" ").join("&nbsp;") + "</span><br>");
            if (outres.evaluate[i].continue) {
                if (outres.evaluate[i].hasOwnProperty("output")) {
                    $("#result_main").append(outres.evaluate[i].output.split("\n").join("<br>") + "<br>");
                }
                createCommandPromptInRconsole("+");
                if (outres.evaluate[i].hasOwnProperty("partial")) {
                    tempcommand += temp[0].substring(outres.evaluate[i].partial.length, temp[0].length) + "\n";
                }
                else {
                    tempcommand += temp.join("\n") + "\n";
                }
            }
            else {
                tempcommand = "";
                if (outres.evaluate[i].hasOwnProperty("message")) {
                    $("#result_main").append("<span style='color:red'>" + outres.evaluate[i].message.split("\n").join("<br>") + "</span>");
                }
                else {
                    $("#result_main").append("<br>");
                }
                if (outres.evaluate[i].hasOwnProperty("output")) {
                    $("#result_main").append(outres.evaluate[i].output.split("\n").join("<br>") + "<br><br>");
                }
                if (outres.evaluate[i].hasOwnProperty("warning")) {
                    $("#result_main").append("<span style='color:red'>" + outres.evaluate[i].warning + "</span><br><br>");
                }
                if (outres.evaluate[i].hasOwnProperty("error")) {
                    $("#result_main").append("<span style='color:red'>" + outres.evaluate[i].error + "</span><br><br>");
                }
                if (outres.evaluate[i].hasOwnProperty("view")) {
                    viewdata = outres.evaluate[i].view;
                }
                if (i == outres.evaluate.length - 1) {
                    createCommandPromptInRconsole();
                }
            }
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            });
            if (focus == "result") {
                $("#tempdiv").click();
            }
        }
        if (viewdata != "") {
            openDataEditor(viewdata);
        }
        reset_outres();
        removeSpinner();
    }
    else {
        if (updatecounter > 12) {
            createSpinner();
        }
        setTimeout(printRcommand, 50);
    }
}
function resizePlot() {
    updatecounter += 1;
    if (updatecounter < 21) {
        if (responseR) {
            updatecounter = 0;
            $("#plotdiv_main").css({
                "background-image": "url('css/images/plot.svg?" + new Date().getTime() + "')", 
                "background-size": "100% 100%"
            });
        }
        else {
            setTimeout(resizePlot, 100);
        }
    }
    else {
        updatecounter = 0;
    }
}
function printWhenOutputChanges() {
    updatecounter += 1;
    if (updatecounter < 21) {
        if (responseR) {
            updatecounter = 0;
            refresh_cols("all");
                if ($("#minimize").length) {
                    checkeqtt();
                }
            history[(histindex < history.length)?history.length:histindex] = string_command;
            histindex = history.length;
            var header = strwrap(string_command, cpl.result, "  ");
            $("#tempdiv").remove();
            $("#result_main").append("<span style='color:#932192'>> </span><span style='color:blue'>" + header + "</span><br>");
            if (objname != "") {
                $("#result_main").append("<span style='color:#932192'>" + "> " + "</span>");
                $("#result_main").append("<span style='color:blue'>" + objname + "</span><br>");
            }
            var toprint = "";
            if (outres.error !== null) {
                for (var i = 0; i < outres.error.length; i++) {
                    toprint += "<br>" + outres.error[i];
                }
                $("#result_main").append("<span style='color:red'>" +
                                         (toprint.split(" ").join("&nbsp;") + "<br>") +
                                         "</span>");
            }
            else {
                if (outres.console !== null) {
                    for (var i = 0; i < outres.console.length; i++) {
                        toprint += "<br>" + outres.console[i];
                    }
                    $("#result_main").append(toprint.split(" ").join("&nbsp;") + "<br>");
                }
                else {
                    $("#result_main").append("<br>");
                }
                if (outres.warning !== null) {
                    toprint = "";
                    for (var i = 0; i < outres.warning.length; i++) {
                        toprint += "<br>" + outres.warning[i];
                    }
                    $("#result_main").append("<span style='color:red'>" +
                                             (toprint.split(" ").join("&nbsp;") + "<br>") +
                                             "</span>");
                }
            }
            createCommandPromptInRconsole();
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1000);
        }
        else {
            setTimeout(printWhenOutputChanges, 50);
        }
    }
    else {
        $("#result_main").append("<br><span style='color:red'>Error in printWhenOutputChanges:<br>R takes too long to respond.</span><br>");
        $("#result_main").animate({
            scrollTop: $("#result_main")[0].scrollHeight
        }, 1000);
        createCommandPromptInRconsole();
        updatecounter = 0;
    }
}
function print_dirs() {
    if ($("#import").length) {
        if (dirfile.filename == "error!") {
            if (papers["import"]["main"].stdir_text.attr("text") != "") {
                papers["import"]["main"].glow.show();
            }
        }
        else {
            papers["import"]["main"].glow.hide();
            papers["import"]["main"].stdir_text.attr({"text": ""});
            papers["import"]["path"].goToDir = function(dir) {
                dirfile_chosen[0] = "dir";
                dirfile_chosen[1] = dir;
                dirfile_chosen[2] = "";
                pathcopy = dirfile.filepath;
                Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                printDirsWhenPathChanges();
            }
            setPath(papers["import"]["path"], dirfile.wd, settings["import"]["inside"]["path"]["width"]);
            var rw = 400; 
            papers["import"]["dirs"].clear();
            var i, aaa, bbb, ccc, toprint, printplus;
            var row = 10;
            var fill_opacity = 0.3;
            var rects_back = papers["import"]["dirs"].set();
            var texts = papers["import"]["dirs"].set();
            var rects = papers["import"]["dirs"].set();
            var pluses = papers["import"]["dirs"].set();
            var x = 30; 
            var dirs_length = 1*((dirfile.dirs === null)?0:dirfile.dirs.length);
            var files_length = 1*((dirfile.files === null)?0:dirfile.files.length);
            var clicked = -1;
            for (i = 0; i < (1 + dirs_length + files_length); i++) {
                printplus = false;
                bbb = papers["import"]["dirs"].rect(0, row - 10, rw, 20);
                bbb.id = i;
                rects_back.push(bbb);
                if (i == 0) {
                    toprint = "..";
                    pluses.push(papers["import"]["dirs"].text(0, 0, ""));
                }
                else {
                    if (i < (dirs_length + 1)) {
                        printplus = true;
                        toprint = dirfile.dirs[i - 1];
                    }
                    else {
                        toprint = dirfile.files[i - dirs_length - 1];
                    }
                }
                if (printplus) {
                    pluses.push(sat(papers["import"]["dirs"].text(x - 20, row, "+")));
                }
                aaa = sat(papers["import"]["dirs"].text(x, row, toprint));
                texts.push(aaa);
                ccc = papers["import"]["dirs"].rect(0, row - 10, rw, 20);
                ccc.id = i;
                ccc.txt = (i == 0)?(".."):(toprint);
                ccc.click(function() {
                        pluses.forEach(function(e) {
                            e.attr({fill: "#000000"});
                        });
                        if (clicked >= 0 && clicked != this.id) {
                            rects_back[clicked].attr({fill: "#e6f2da", "stroke-opacity": 0, "fill-opacity": 1 - 1*(clicked % 2 === 0)});
                            texts[clicked].attr({fill: "#000000"});
                        }
                        clicked = this.id;
                        rects_back[this.id].attr({fill: "#79a74c", "stroke-opacity": 0, "fill-opacity": 1});
                        texts[this.id].attr({fill: "#ffffff"});
                        if (this.id < (dirs_length + 1)) {
                            pluses[this.id].attr({fill: "#ffffff"});
                        }
                        dirfile_chosen[0] = (this.id < (dirs_length + 1))?"dir":"file";
                        dirfile_chosen[1] = (this.txt == "..")?((dirfile_chosen[1] == "..")?"...":".."):this.txt;
                        if (dirfile_chosen[0] == "file") {
                            tempdatainfo.rownames = "foobar";
                            tempdatainfo.nrows = 0;
                            commobj.read_table.counter += 1;
                            Shiny.onInputChange("read_table", commobj.read_table);
                            papers["import"]["cols"].clear();
                            $(papers["import"]["cols"].canvas).width(100);
                            papers["import"]["cols"].text(10, 11, "Loading...").attr({"text-anchor": "start", "font-size": "14px"});
                            updatecounter = 0;
                            checkIfDataLoadedInR();
                            pathcopy = dirfile.filepath;
                            updatecounter2 = 0;
                            Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                            consoleIfPathChanges();
                        }
                    })
                    .dblclick(function(event) {
                        dirfile_chosen[0] = (this.id < (dirs_length + 1))?"dir":"file";
                        dirfile_chosen[1] = (this.txt == "..")?((dirfile_chosen[1] == "..")?"...":".."):this.txt;
                        dirfile_chosen[2] = "";
                        papers["import"]["main"].stdir_text.attr({"text": ""});
                        if (dirfile_chosen[0] == "dir") {
                            pathcopy = dirfile.filepath;
                            Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                            printDirsWhenPathChanges();
                        }
                    });
                rects.push(ccc);
                row += 20;
            }
            rects_back.forEach(function(e) {
                e.attr({fill: (e.id % 2 === 0)?"#ffffff":"#e6f2da", "stroke-opacity": 0, "fill-opacity": 1});
            });
            rects.forEach(function(e) {
                e.attr({fill: "#ffffff", "stroke-opacity": 0, "fill-opacity": 0});
            });
            rects_back.toBack();
            rects.toFront();
            canvas_height = Math.max(400, rects.getBBox().height);
            $(papers["import"]["dirs"].canvas).height(canvas_height);
            $("#import_dirs").css({height: canvas_height});
            if (dirfile_chosen[0] == "dir") {
                $("#import_dirs").scrollTop(0);
            }
        }
    }
    if ($("#export").length) {
        if (dirfile.filename == "error!") {
            if (papers["export"]["main"].stdir_text.attr("text") != "") {
                papers["export"]["main"].glow.show();
            }
        }
        else {
            papers["export"]["main"].glow.hide();
            papers["export"]["main"].stdir_text.attr({"text": ""});
            papers["export"]["path"].goToDir = function(dir) {
                dirfile_chosen[0] = "dir";
                dirfile_chosen[1] = dir;
                pathcopy = dirfile.filepath;
                Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                printDirsWhenPathChanges();
            }
            setPath(papers["export"]["path"], dirfile.wd, settings["export"]["inside"]["path"]["width"]);
            var rw = 400; 
            papers["export"]["dirs"].clear();
            var i, aaa, bbb, ccc, extoprint, printplus;
            var row = 10;
            var fill_opacity = 0.3;
            var exrects_back = papers["export"]["dirs"].set();
            var extexts = papers["export"]["dirs"].set();
            var exrects = papers["export"]["dirs"].set();
            var expluses = papers["export"]["dirs"].set();
            var x = 30; 
            var exdirs_length = 1*((dirfile.dirs === null)?0:dirfile.dirs.length);
            var exfiles_length = 1*((dirfile.files === null)?0:dirfile.files.length);
            var exclicked = -1;
            papers["export"]["main"].ovr.hideIt();
            if (dirfile.files.length != 0) {
                if (dirfile.files.indexOf(commobj["export"].filename) >= 0) {
                    papers["export"]["main"].ovr.showIt();
                }
            }
            console_command("export");
            for (i = 0; i < (1 + exdirs_length + exfiles_length); i++) {
                printplus = false;
                bbb = papers["export"]["dirs"].rect(0, row - 10, rw, 20);
                bbb.id = i;
                exrects_back.push(bbb);
                if (i == 0) {
                    extoprint = "..";
                    expluses.push(papers["export"]["dirs"].text(0, 0, ""));
                }
                else {
                    if (i < (exdirs_length + 1)) {
                        printplus = true;
                        extoprint = dirfile.dirs[i - 1];
                    }
                    else {
                        extoprint = dirfile.files[i - exdirs_length - 1];
                    }
                }
                if (printplus) {
                    expluses.push(papers["export"]["dirs"].text(x - 20, row, "+").attr({"text-anchor": "start", "font-size": "14px"}));
                }
                aaa = sat(papers["export"]["dirs"].text(x, row, extoprint));
                extexts.push(aaa);
                ccc = papers["export"]["dirs"].rect(0, row - 10, rw, 20);
                ccc.id = i;
                ccc.txt = (i == 0)?(".."):(extoprint);
                ccc.click(function() {
                        expluses.forEach(function(e) {
                            e.attr({fill: "#000000"});
                        });
                        if (exclicked >= 0 && exclicked != this.id) {
                            exrects_back[exclicked].attr({fill: "#e6f2da", "stroke-opacity": 0, "fill-opacity": 1*(exclicked % 2 === 0)});
                            extexts[exclicked].attr({fill: "#000000"});
                        }
                        exclicked = this.id;
                        exrects_back[this.id].attr({fill: "#79a74c", "stroke-opacity": 0, "fill-opacity": 1});
                        extexts[this.id].attr({fill: "#ffffff"});
                        if (this.id < (exdirs_length + 1)) {
                            expluses[this.id].attr({fill: "#ffffff"});
                        }
                        if (this.id > exdirs_length) { 
                            papers["export"]["main"].newname.attr({"text": this.txt});
                            commobj["export"].filename = this.txt;
                            papers["export"]["main"].ovr.showIt();
                        }
                    })
                    .dblclick(function() {
                        dirfile_chosen[0] = (this.id < (exdirs_length + 1))?"dir":"file";
                        dirfile_chosen[1] = (this.txt == "..")?((dirfile_chosen[1] == "..")?"...":".."):this.txt;
                        if (dirfile_chosen[0] == "dir") {
                            pathcopy = dirfile.filepath;
                            Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                            consoleIfPathChanges();
                            printDirsWhenPathChanges();
                        }
                    });
                exrects.push(ccc);
                row += 20;
            }
            exrects_back.forEach(function(e) {
                e.attr({fill: (e.id % 2 === 0)?"#e6f2da":"#ffffff", "stroke-opacity": 0, "fill-opacity": 1});
            });
            exrects.forEach(function(e) {
                e.attr({fill: "#ffffff", "stroke-opacity": 0, "fill-opacity": 0});
            });
            exrects_back.toBack();
            exrects.toFront();
            canvas_height = Math.max(400, exrects.getBBox().height);
            $(papers["export"]["dirs"].canvas).height(canvas_height);
            $("#export_dirs").css({height: canvas_height});
                $("#export_dirs").scrollTop(0);
        }
    }
    if ($("#saveRplot").length) {
        if (dirfile.filename == "error!") {
            if (papers["saveRplot"]["main"].stdir_text.attr("text") != "") {
                papers["saveRplot"]["main"].glow.show();
            }
        }
        else {
            papers["saveRplot"]["main"].glow.hide();
            papers["saveRplot"]["main"].stdir_text.attr({"text": ""});
            papers["saveRplot"]["path"].goToDir = function(dir) {
                dirfile_chosen[0] = "dir";
                dirfile_chosen[1] = dir;
                pathcopy = dirfile.filepath;
                Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                printDirsWhenPathChanges();
            }
            setPath(papers["saveRplot"]["path"], dirfile.wd, settings["saveRplot"]["inside"]["path"]["width"]);
            var rw = 400; 
            papers["saveRplot"]["dirs"].clear();
            if (savetexts !== void 0) {
                savetexts.remove();
            }
            var aaa, bbb, ccc, svtoprint, printplus;
            var row = 10;
            var fill_opacity = 0.3;
            var saverects_back = papers["saveRplot"]["dirs"].set();
            var savetexts = papers["saveRplot"]["dirs"].set();
            var saverects = papers["saveRplot"]["dirs"].set();
            var savepluses = papers["saveRplot"]["dirs"].set();
            var x = 30; 
            var savedirs_length = 1*((dirfile.dirs === null)?0:dirfile.dirs.length);
            var savefiles_length = 1*((dirfile.files === null)?0:dirfile.files.length);
            var svclicked = -1;
            papers["saveRplot"]["main"].ovr.hideIt();
            if (dirfile.files.length != 0) {
                if (dirfile.files.indexOf((commobj["saveRplot"].filename + "." + commobj["saveRplot"].type)) >= 0) {
                    papers["saveRplot"]["main"].ovr.showIt();
                }
            }
            console_command("saveRplot");
            for (var i = 0; i < (1 + savedirs_length + savefiles_length); i++) {
                printplus = false;
                bbb = papers["saveRplot"]["dirs"].rect(0, row - 10, rw, 20);
                bbb.id = i;
                saverects_back.push(bbb);
                if (i == 0) {
                    svtoprint = "..";
                    savepluses.push(papers["saveRplot"]["dirs"].text(0, 0, ""));
                }
                else {
                    if (i < (savedirs_length + 1)) {
                        printplus = true;
                        svtoprint = dirfile.dirs[i - 1];
                    }
                    else {
                        svtoprint = dirfile.files[i - savedirs_length - 1];
                    }
                }
                if (printplus) {
                    savepluses.push(papers["saveRplot"]["dirs"].text(x - 20, row, "+").attr({"text-anchor": "start", "font-size": "14px"}));
                }
                aaa = sat(papers["saveRplot"]["dirs"].text(x, row, svtoprint));
                savetexts.push(aaa);
                ccc = papers["saveRplot"]["dirs"].rect(0, row - 10, rw, 20);
                ccc.id = i;
                ccc.txt = (i == 0)?(".."):(svtoprint);
                ccc.click(function() {
                        savepluses.forEach(function(e) {
                            e.attr({fill: "#000000"});
                        });
                        if (svclicked >= 0 && svclicked != this.id) {
                            saverects_back[svclicked].attr({fill: "#e6f2da", "stroke-opacity": 0, "fill-opacity": 1*(svclicked % 2 === 0)});
                            savetexts[svclicked].attr({fill: "#000000"});
                        }
                        svclicked = this.id;
                        saverects_back[this.id].attr({fill: "#79a74c", "stroke-opacity": 0, "fill-opacity": 1});
                        savetexts[this.id].attr({fill: "#ffffff"});
                        if (this.id < (savedirs_length + 1)) {
                            savepluses[this.id].attr({fill: "#ffffff"});
                        }
                        if (this.id > savedirs_length) { 
                            papers["saveRplot"]["main"].filename.attr({"text": this.txt});
                            commobj.saveRplot.filename = this.txt;
                            papers["saveRplot"]["main"].ovr.showIt();
                        }
                    })
                    .dblclick(function() {
                        dirfile_chosen[0] = (this.id < (savedirs_length + 1))?"dir":"file";
                        dirfile_chosen[1] = (this.txt == "..")?((dirfile_chosen[1] == "..")?"...":".."):this.txt;
                        if (dirfile_chosen[0] == "dir") {
                            pathcopy = dirfile.filepath;
                            Shiny.onInputChange("dirfile_chosen", dirfile_chosen);
                            consoleIfPathChanges();
                            printDirsWhenPathChanges();
                        }
                    });
                saverects.push(ccc);
                row += 20;
            }
            saverects_back.forEach(function(e) {
                e.attr({fill: (e.id % 2 === 0)?"#e6f2da":"#ffffff", "stroke-opacity": 0, "fill-opacity": 1});
            });
            saverects.forEach(function(e) {
                e.attr({fill: "#ffffff", "stroke-opacity": 0, "fill-opacity": 0});
            });
            saverects_back.toBack();
            saverects.toFront();
            canvas_height = Math.max(400, saverects.getBBox().height);
            $(papers["saveRplot"]["dirs"].canvas).height(canvas_height);
            $("#saveRplot_dirs").css({height: canvas_height});
                $("#saveRplot_dirs").scrollTop(0);
        }
    }
}
function updateWhenDataChanged() {
    updatecounter += 1;
    if (updatecounter < 21) { 
        if (responseR) {
            print_data();
            updatecounter = 0; 
        }
        else {
            setTimeout(updateWhenDataChanged, 50);
        }
    }
    else {
        updatecounter = 0; 
    }
}
function printIfDirsFilesChange() {
    updatecounter += 1;
    if (updatecounter < 21) {
        var test = "";
        if (dirfile.dirs != null) {
            test += dirfile.dirs.toString();
        }
        if (dirfile.files != null) {
            test += dirfile.files.toString();
        }
        if (test != dirsfilescopy) {
            print_dirs();
            updatecounter = 0;
        }
        else {
            setTimeout(printIfDirsFilesChange, 50);
        }
    }
    else {
        updatecounter = 0;
    }
}
function printDirsWhenPathChanges() {
    updatecounter += 1;
    if (updatecounter < 21) {
        if (dirfile.filepath != pathcopy) {
            print_dirs();
        }
        else {
            setTimeout(printDirsWhenPathChanges, 50);
        }
    }
    else {
        updatecounter = 0;
    }
}
function textbox(paper, settings) {
    var x = settings.x;
    var y = settings.y;
    var text = settings.text;
    var width = settings.width;
    var height = settings.height;
    var result = {};
    result.text = sat(paper.text(x, y, text), {"clip": (x - 5) + ", " + (y - 10) + ", " + width + ", " + height});
    result.rect = sat(paper.rect(x - 5, y - 10, width, height, 3));
    return(result);
}
var lastX, absoluteX, newpos;
function dragStart(slider) {
    return function() {
        lastX = 0;
        var getBB = this.getBBox()
        absoluteX = getBB.x + getBB.width/2;
        if (this.id == "thsetter" && commobj.calibrate.type == "crisp") {
            commobj.calibrate.findth = false;
            papers["calibrate"]["main"]["findth"].uncheck();
        }
    }
};
function dragMove(slider) {
    return function(dx, dy) {
        var newX = dx - lastX;
        if (absoluteX + dx > this.right) {
            newX = this.right - absoluteX - lastX;
        }
        if (absoluteX + dx < this.left) {
            newX = this.left - absoluteX - lastX;
        }
        this.translate(newX, 0);
        if (absoluteX + dx < this.left) {
            lastX = this.left - absoluteX;
        }
        else {
            lastX += newX;
        }
        newpos = (absoluteX + lastX - this.left)/(this.right - this.left);
        newpos = this.min + newpos*(this.max - this.min);
        if (this.id == "thsetter") {
            newpos = Math.round(newpos*1000)/1000;
            ths[this.name].attr({"text": newpos});
            commobj.calibrate.thresholds[this.name] = newpos;
        }
        else if (this.id == "xyplot") {
            papers["xyplot"]["main"].labelRotation = newpos;
            if (xyplotdata.length > 0) {
                createLabels(papers["xyplot"]["main"]);
            }
        }
    }
};
function dragStop(slider) {
    return function() {
        if (this.id == "thsetter") {
            var paper = papers["calibrate"]["main"];
            var crisp = commobj.calibrate.type == "crisp";
            var bell = !commobj.calibrate.end;
            var sort = crisp;
            if (commobj.calibrate.type == "fuzzy") {
                var filled = 0;
                for (var j = 0; j < 6; j++) {
                    if (commobj.calibrate.thresholds[j] != "") {
                        filled += 1;
                    }
                }
                sort = sort || (bell)?(filled == 6):(filled == 3);
            }
            if (sort) {
                commobj.calibrate.thresholds = sortArray(commobj.calibrate.thresholds);
                for (var j = 0; j < 6; j++) {
                    ths[j].attr({"text": commobj.calibrate.thresholds[j]});
                }
            }
            if (crisp) {
                commobj.calibrate.thscopycrp = copy(commobj.calibrate.thresholds);
            }
            else { 
                if (bell) {
                    commobj.calibrate.thscopyfuzb = copy(commobj.calibrate.thresholds);
                }
                else {
                    commobj.calibrate.thscopyfuzs = copy(commobj.calibrate.thresholds);
                }
            }
            if (commobj.calibrate.type == "crisp" || (commobj.calibrate.type == "fuzzy" && !checkfuzzy())) {
                drawPointsAndThresholds();
            }
            console_command("calibrate");
        }
        else if (this.id == "xyplot") {
            papers["xyplot"]["main"].labelRotation = newpos;
        }
        if (colclicks.calibrate !== void 0) {
            commobj.calibrate.x = getTrueKeys(colclicks.calibrate.x)[0];
            if (commobj.calibrate.x === void 0) {
                commobj.calibrate.x = "";
            }
        }
        if (commobj.calibrate.x != "") {
            if (commobj.calibrate.type == "fuzzy" && checkfuzzy()) {
                commobj.calibrate.counter += 1;
                commobj.calibrate.thsetter = true;
                Shiny.onInputChange("calibrate", commobj.calibrate);
            }
        }
    }
};
var lastY, absoluteY;
function dragSortStart(sortoption) {
    return function() {
        lastY = 0;
        var getBB = sortoption[2].getBBox();
        absoluteY = getBB.y + getBB.height/2;
        sortoption.toFront();
    }
};
function dragSortMove(sortoption) {
    return function(dx, dy) {
        var newY = dy - lastY;
        if (absoluteY + dy > sortoption[2].bottom) {
            newY = sortoption[2].bottom - absoluteY - lastY;
        }
        if (absoluteY + dy < sortoption[2].top) {
            newY = sortoption[2].top - absoluteY - lastY;
        }
        sortoption.translate(0, newY);
        if (absoluteY + dy < sortoption[2].top) {
            lastY = sortoption[2].top - absoluteY;
        }
        else {
            lastY += newY;
        }
    }
};
function dragSortStop(sortoption) {
    return function() {
        var newBB = sortoption[2].getBBox();
        var middle = newBB.y + newBB.height/2;
        var oldposition, newposition, decid;
        for (var i = 0; i < 3; i++) {
            if (absoluteY > papers["tt"]["main"].coordsy[i]) {
                oldposition = i;
            }
            if (middle > papers["tt"]["main"].coordsy[i]) {
                newposition = i;
            }
        }
        if (oldposition == newposition) {
            if (commobj.tt.sort_sel[sortoption[2].name]) {
                sortoption[0].attr({fill: "#eeeeee", stroke: "none"});
                sortoption[1].attr({fill: "black", "text-anchor": "start", "font-size": "14px"});
            }
            else {
                sortoption[0].attr({fill: "#79a74c", stroke: "none"});
                sortoption[1].attr({fill: "white", "text-anchor": "start", "font-size": "14px"});
            }
            commobj.tt.sort_sel[sortoption[2].name] = !commobj.tt.sort_sel[sortoption[2].name];
            sortoption[0].backcolor = commobj.tt.sort_sel[sortoption[2].name];
            sortoption.translate(0, papers["tt"]["main"].coordsy[oldposition] - newBB.y);
        }
        else {
            var positions = copy(papers["tt"]["main"].positions);
            var distomove;
            distomove = papers["tt"]["main"].coordsy[oldposition] - papers["tt"]["main"].sortsets[positions[newposition]].getBBox().y;
            papers["tt"]["main"].sortsets[positions[newposition]].translate(0, distomove);
            sortoption.translate(0, papers["tt"]["main"].coordsy[newposition] - newBB.y);
            papers["tt"]["main"].positions[oldposition] = positions[newposition];
            papers["tt"]["main"].positions[newposition] = positions[oldposition];
            if (Math.abs(newposition - oldposition) == 2) {
                distomove = papers["tt"]["main"].coordsy[oldposition] - papers["tt"]["main"].sortsets[positions[1]].getBBox().y;
                papers["tt"]["main"].sortsets[positions[1]].translate(0, distomove);
                distomove = papers["tt"]["main"].coordsy[1] - papers["tt"]["main"].sortsets[positions[newposition]].getBBox().y;
                papers["tt"]["main"].sortsets[positions[newposition]].translate(0, distomove);
                papers["tt"]["main"].positions[oldposition] = positions[1];
                papers["tt"]["main"].positions[1] = positions[newposition];
            }
            commobj.tt.sort_by = reorder(commobj.tt.sort_by, oldposition, newposition);
            commobj.tt.sort_sel = reorder(commobj.tt.sort_sel, oldposition, newposition);
        }
        var keys = getKeys(commobj.tt.sort_by);
        for (var i = 0; i < 3; i++) {
            papers["tt"]["main"].decrease[i].cover.name = keys[i];
            if (commobj.tt.sort_by[keys[i]]) {
                papers["tt"]["main"].decrease[i].check();
            }
            else {
                papers["tt"]["main"].decrease[i].uncheck();
            }
            if (commobj.tt.sort_sel[keys[i]]) {
                papers["tt"]["main"].decrease[i].showIt();
            }
            else {
                papers["tt"]["main"].decrease[i].hideIt();
            }
        }
        if (getTrueKeys(commobj.tt.sort_sel).length == 0) {
            papers["tt"]["main"].decr.hide();
        }
        else {
            papers["tt"]["main"].decr.show();
        }
        console_command("tt");
    }
}
$("#menu_load").click(function() {
    pingit();
    if ($("#load").length) {
        showDialogToFront(settings["load"]);
    }
    else {
        current_command = "load";
        createDialog(settings["load"]);
        makePapers(settings["load"], papers);
        draw_load(papers["load"]["main"]);
        refresh_cols("load");
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_import").click(function() {
    pingit();
    if ($("#import").length) {
        showDialogToFront(settings["import"]);
    }
    else {
        current_command = "import";
        createDialog(settings["import"]);
        makePapers(settings["import"], papers);
        $(papers["import"]["cols"].canvas).height(20);
        draw_import(papers["import"]["main"]);
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_export").click(function() {
    pingit();
    if ($("#export").length) {
        showDialogToFront(settings["export"]);
    }
    else {
        current_command = "export";
        createDialog(settings["export"]);
        makePapers(settings["export"], papers);
        refresh_cols("export");
        draw_export(papers["export"]["main"]);
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
function openDataEditor(dataset) {
    pingit();
    if ($("#data_editor").length) {
        showDialogToFront(settings["data_editor"]);
        click_col("data_editor", "dataset", dataset, others = false);
        print_data();
    }
    else {
        var vscrollbar = false, hscrollbar = false;
        if (info["data"] !== null && dataset != "") {
            vscrollbar = info["data"][dataset].nrows - 1 > visiblerows;
            hscrollbar = info["data"][dataset].ncols - 1 > visiblecols;
        }
        createDialog(settings["data_editor"]);
        makePapers(settings["data_editor"], papers);
        $("#data_editor").width((visiblecols + 1 + 1)*70 + 1*(vscrollbar?scrollbarsWH:0));
        $("#data_editor").height((visiblerows + 1 + 2)*20 + 1*(hscrollbar?scrollbarsWH:0));
        $("#data_editor").draggable({
            drag: function(ev, ui) {
                testX = ui.position.left;
                testY = ui.position.top;
            }
        });
        if (info["data"] !== null) {
            if (dataset != "") {
                commobj["data_editor"].dataset = dataset;
                click_col("data_editor", "dataset", dataset, others = false);
            }
            scrollobj.scrollvh = {};
            var datasets = getKeys(info["data"]);
            for (var i = 0; i < datasets.length; i++) {
                scrollobj.scrollvh[datasets[i]] = info["data"][datasets[i]].scrollvh;
            }
            scrollobj.dataset = dataset;
            scrollobj.counter += 1;
            scrollobj.visiblerows = visiblerows;
            scrollobj.visiblecols = visiblecols;
            scrollobj.alldata = true;
            responseR = false;
            Shiny.onInputChange("scrollobj", scrollobj);
            updatecounter = 0;
            responseR = false;
            updateWhenDataChanged();
        }
        else {
            print_data();  
        }
        $("#data_editor").resizable({
            start: function(event, ui) {
            },
            resize: function(event, ui) {
                if (info["data"] !== null) {
                    vscrollbar = info["data"][commobj["data_editor"].dataset].nrows - 1 > visiblerows;
                    hscrollbar = info["data"][commobj["data_editor"].dataset].ncols - 1 > visiblecols;
                }
                $("#data_editor_main").width($("#data_editor").width());
                $("#data_editor_main").height($("#data_editor").height() - 20);
                var width = $("#data_editor").width() - 70;
                var height = $("#data_editor").height() - 40;
                $("#data_editor_body").width(width);
                $("#data_editor_body").height(height);
                $("#data_editor_colnames").width(width - 1*(vscrollbar?scrollbarsWH:0));
                $("#data_editor_rownames").height(height - 1*(hscrollbar?scrollbarsWH:0));
                visiblerows = Math.round($("#data_editor_rownames").height()/20) - 1; 
                visiblecols = Math.round($("#data_editor_colnames").width() /70) - 1;
                $(papers["data_editor"]["body"].canvas).width(width - 15);
                $(papers["data_editor"]["body"].canvas).height(height - 15);
                $(papers["data_editor"]["colnames"].canvas).width(width - 1*(vscrollbar?scrollbarsWH:0) - 15);
                $(papers["data_editor"]["rownames"].canvas).height(height - 1*(vscrollbar?scrollbarsWH:0) - 15);
            },
            stop: function(event, ui) {
                var dataset = commobj["data_editor"].dataset;
                if (info["data"] !== null) {
                    vscrollbar = info["data"][dataset].nrows - 1 > visiblerows;
                    hscrollbar = info["data"][dataset].ncols - 1 > visiblecols;
                }
                $("#data_editor").width((visiblecols + 2)*70 + 1*(vscrollbar?scrollbarsWH:0));
                $("#data_editor").height((visiblerows + 3)*20 + 1*(hscrollbar?scrollbarsWH:0));
                $("#data_editor_main").width($("#data_editor").width());
                $("#data_editor_main").height($("#data_editor").height() - 20);
                if (dataset != "") {
                    var svh = info["data"][dataset].scrollvh.toString();
                    info["data"][dataset].scrollvh[2] = info["data"][dataset].scrollvh[0] + visiblerows;
                    info["data"][dataset].scrollvh[3] = info["data"][dataset].scrollvh[1] + visiblecols;
                    if (info["data"][dataset].scrollvh.toString() !== svh && info["data"] !== null) {
                        visibledata = info["data"][dataset].theData.toString();
                        coordscopy = info["data"][dataset].dataCoords;
                        scrollobj.scrollvh = {};
                        var datasets = getKeys(info["data"]);
                        for (var i = 0; i < datasets.length; i++) {
                            scrollobj.scrollvh[datasets[i]] = info["data"][datasets[i]].scrollvh;
                        }
                        scrollobj.dataset = dataset;
                        scrollobj.counter += 1;
                        scrollobj.visiblerows = visiblerows;
                        scrollobj.visiblecols = visiblecols;
                        scrollobj.alldata = true;
                        responseR = false;
                        Shiny.onInputChange("scrollobj", scrollobj);
                        updatecounter = 0;
                        updateWhenDataChanged();
                    }
                    else {
                        print_data();
                    }
                }
                else {
                    print_data();
                    var width = $("#data_editor").width() - 70;
                    var height = $("#data_editor").height() - 40;
                    $(papers["data_editor"]["body"].canvas).width(width);
                    $(papers["data_editor"]["body"].canvas).height(height);
                    $(papers["data_editor"]["colnames"].canvas).width(width - 1*(vscrollbar?scrollbarsWH:0));
                    $(papers["data_editor"]["rownames"].canvas).height(height - 1*(vscrollbar?scrollbarsWH:0));
                }
            }
        });
        $("#data_editor_body").scroll(function () {
            $("#data_editor_rownames").scrollTop($("#data_editor_body").scrollTop());
            $("#data_editor_colnames").scrollLeft($("#data_editor_body").scrollLeft());
            clearTimeout($.data(this, 'scrollCheck'));
            $.data(this, "scrollCheck", setTimeout(function() {
                var vertical = $("#data_editor_body").scrollTop();
                var horizontal = $("#data_editor_body").scrollLeft(); 
                var cellstoright = Math.round(horizontal/70);
                var cellsdown = Math.round(vertical/20);
                var dataset = commobj["data_editor"].dataset;
                scrolleftop[dataset] = {"vertical": vertical, "horizontal": horizontal}
                var change = false;
                if (cellstoright != info["data"][dataset].scrollvh[1] && info["data"] !== null) {
                    change = true;
                    info["data"][dataset].scrollvh[1] = cellstoright;
                }
                if (cellsdown != info["data"][dataset].scrollvh[0] && info["data"] !== null) {
                    change = true;
                    info["data"][dataset].scrollvh[0] = cellsdown;
                }
                if (change && info["data"] !== null) {
                    visibledata = info["data"][dataset].theData.toString();
                    coordscopy = info["data"][dataset].dataCoords;
                    scrollobj.scrollvh = {};
                    scrollobj.scrollvh[dataset] = info["data"][dataset].scrollvh;
                    scrollobj.dataset = dataset;
                    scrollobj.counter += 1;
                    scrollobj.visiblerows = visiblerows;
                    scrollobj.visiblecols = visiblecols;
                    scrollobj.alldata = false;
                    responseR = false;
                    Shiny.onInputChange("scrollobj", scrollobj);
                    updateWhenDataChanged();
                }
            }, 100));
        });
    }
}
$("#menu_data_editor").click(function() {
    openDataEditor(commobj["data_editor"].dataset);
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_calibrate").click(function() {
    pingit();
    if ($("#calibrate").length) {
        showDialogToFront(settings["calibrate"]);
    }
    else {
        createDialog(settings["calibrate"]);
        makePapers(settings["calibrate"], papers);
        thsetter_content = papers["calibrate"]["main"].set();
        refresh_cols("calibrate");
        draw_calib(papers["calibrate"]["main"]);
        $("#calibrate").resizable({
            start: function() {
                showDialogToFront(settings["calibrate"]);
            },
            resize: function() {
                $(this).height(settings["calibrate"].height);
                $("#calibrate_main").width($("#calibrate").width());
                $(papers["calibrate"]["main"].canvas).width($("#calibrate").width());
                papers["calibrate"]["main"].thsetter_frame.attr({width: $("#calibrate").width() - 175});
                papers["calibrate"]["main"].Run.transform("t" + ($("#calibrate").width() - 500) + ",0");
            },
            stop: function () {
                if (poinths.vals.length > 0) {
                    drawPointsAndThresholds();
                }
            }
        });
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_recode").click(function() {
    pingit();
    if ($("#recode").length) {
        showDialogToFront(settings["recode"]);
    }
    else {
        createDialog(settings["recode"]);
        makePapers(settings["recode"], papers);
        refresh_cols("recode");
        draw_recode(papers["recode"]["main"]);
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_create_tt").click(function() {
    pingit();
    if ($("#tt").length) {
        showDialogToFront(settings["tt"]);
    }
    else {
        createDialog(settings["tt"]);
        makePapers(settings["tt"], papers);
        refresh_cols("tt");
        draw_tt(papers["tt"]["main"]);
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_find_rows").click(function() {
    pingit();
    if ($("#findRows").length) {
        showDialogToFront(settings["findRows"]);
    }
    else {
        createDialog(settings["findRows"]);
        makePapers(settings["findRows"], papers);
        refresh_cols("findRows");
        draw_findRows(papers["findRows"]["main"]);
        current_command = "findRows";
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_minimize").click(function() {
    pingit();
    if ($("#minimize").length) {
        showDialogToFront(settings["minimize"]);
    }
    else {
        createDialog(settings["minimize"]);
        makePapers(settings["minimize"], papers);
        draw_minimize(papers["minimize"]["main"]);
        refresh_cols("minimize");
        checkeqtt()
        current_command = "minimize";
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_xyplot").click(function() {
    pingit();
    if ($("#xyplot").length) {
        showDialogToFront(settings["xyplot"]);
    }
    else {
        createDialog(settings["xyplot"]);
        makePapers(settings["xyplot"], papers);
        refresh_cols("xyplot");
        draw_xyplot(papers["xyplot"]["main"]);
        $("#xyplot").resizable({
            resize: function () {
                var paper = papers["xyplot"]["main"];
                $("#xyplot").width($("#xyplot_main").height() + 173);
                $(paper.canvas).width($("#xyplot").width());
                $(paper.canvas).height($("#xyplot").height() - 20);
                $("#xyplot_main").width($("#xyplot").width());
                $("#xyplot_main").height($("#xyplot").height() - 20);
            },
            stop: function() {
                var paper = papers["xyplot"]["main"];
                $(paper.canvas).width($("#xyplot").width() - 50);
                $(paper.canvas).height($("#xyplot").height() - 70);
                $(paper.canvas).width($("#xyplot").width());
                $(paper.canvas).height($("#xyplot").height() - 20);
                paper.scale = Math.min(($(paper.canvas).width() - paper.sx - 10)/paper.dim, ($(paper.canvas).height() - paper.sy - 47)/paper.dim);
                draw_xyplot(paper);
            }
        });
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_venn").click(function() {
    pingit();
    if ($("#venn").length) {
        showDialogToFront(settings["venn"]);
    }
    else {
        createDialog(settings["venn"]);
        makePapers(settings["venn"], papers);
        draw_venn(papers["venn"]["main"]);
        $("#venn").resizable({
            resize: function() {
                if ($(this).height() < settings["venn"].minHeight) {
                    $(this).height(settings["venn"].minHeight);
                }
                if ($(this).width() < settings["venn"].minWidth) {
                    $(this).width(settings["venn"].minWidth);
                }
                var paper = papers["venn"]["main"];
                paper.hover = false;
                $(this).height($(this).width() + 70);
                $("#venn_main").width($(this).width());
                $("#venn_main").height($(this).height() - 20);
                $(paper.canvas).width($(this).width());
                $(paper.canvas).height($(this).height() - 20);
            },
            stop: function() {
                var paper = papers["venn"]["main"];
                $(paper.canvas).width($(this).width() - 50);
                $(paper.canvas).height($(this).height() - 70);
                $(paper.canvas).width($(this).width());
                $(paper.canvas).height($(this).height() - 20);
                paper.scale = (Math.min($(paper.canvas).width() - 20, $(paper.canvas).height() - 70))/1000;
                draw_venn(paper);
                paper.hover = true;
            }
        });
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_saveRplot").click(function() {
    pingit();
    if ($("#saveRplot").length) {
        showDialogToFront(settings["saveRplot"]);
    }
    else {
        current_command = "saveRplot";
        createDialog(settings["saveRplot"]);
        makePapers(settings["saveRplot"], papers);
        draw_saveRplot(papers["saveRplot"]["main"]);
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_about").click(function() {
    pingit();
    if ($("#about").length) {
        showDialogToFront(settings["about"]);
    }
    else {
        createDialog(settings["about"]);
        var messages = [
            "R package: QCA, version 3.22.5",
            "",
            "Author: Adrian Dușa (dusa.adrian@unibuc.ro)",
            "Former coauthors:",
            "        Alrik Thiem (package QCA version 1.0-0 up to version 1.1-3)",
            "",
            "Contributors:",
            "            jQuery Foundation  (jQuery library and jQuery UI library)",
            "            jQuery contributors (jQuery library and jQuery UI library)",
            "            Vasil Dinkov (smartmenus.js library)",
            "            Dmitry Baranovskiy (raphael.js library)",
            "            Emmanuel Quentin (raphael.inline_text_editing.js library)",
            "            Jimmy Breck-McKye (raphael-paragraph.js library)",
            "",
            "An extensive set of functions to perform Qualitative Comparative",
            "Analysis: crisp sets ('csQCA'), temporal ('tQCA'), multivalue ('mvQCA')",
            "and fuzzy sets ('fsQCA'), using a GUI - graphical user interface.",
            "'QCA' is a methodology that bridges the qualitative and quantitative",
            "divide in social science research. It uses a Boolean algorithm that",
            "results in a minimal causal combination which explains a given",
            "phenomenon.",
            "",
            "To cite this package in publications, please use:",
            "        Dușa, Adrian (2019) QCA with R. A Comprehensive Resource.",
            "        Springer International Publishing."
        ];
        for (var i = 0; i < messages.length; i++) {
            var text = strwrap(messages[i], 80, "  ");
            $("#about_main").append(text + "<br>");
        }
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_changes").click(function() {
    pingit();
    if ($("#changes").length) {
        showDialogToFront(settings["changes"]);
    }
    else {
        createDialog(settings["changes"]);
        function getChanges() {
            updatecounter += 1;
            if (updatecounter < 21) {
                if (responseR) {
                    updatecounter = 0;
                    for (var i = 0; i < outres.length; i++) {
                        var text = strwrap(outres[i].replace("      ", "       "), 87, "  ");
                        $("#changes_main").append(text + "<br>");
                    }
                }
                else {
                    setTimeout(getChanges, 50);
                }
            }
        }
        changes += 1;
        updatecounter = 0;
        outres = new Array();
        responseR = false;
        Shiny.onInputChange("changes", changes);
        getChanges();
    }
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_help").click(function() {
    pingit();
    help += 1;
    Shiny.onInputChange("help", help);
    $("#main_menu").smartmenus('menuHideAll');
    return false;
});
$("#menu_quit").click(function() {
    pingit();
    var quit = 0;
    quit += 1;
    if (confirm("Close Window?")) {
        Shiny.onInputChange("quit", quit);
    }
});
createDialog(settings["command"]);
createDialog(settings["result"]);
createDialog(settings["plotdiv"]);
$("#command").resizable({
    stop: function() {
        var width = $("#command").width() - 14;
        cpl.command = Math.floor(width*1.64/14);
        $("#command_main").html(strwrap(string_command, cpl.command));
    }
});
$("#result").resizable({
    stop: function() {
        var width = $("#result").width() - 25;
        cpl.result = Math.floor(width*1.64/14);
    }
});
$("#plotdiv").resizable({
    stop: function() {
        Rcommand.plotsize[0] = $('#plotdiv_main').width()/96
        Rcommand.plotsize[1] = $('#plotdiv_main').height()/96;
        outres = new Array();
        responseR = false;
        Shiny.onInputChange("plotsize", Rcommand.plotsize);
        updatecounter = 0;
        resizePlot();
    }
});
$("#plotdiv").hide();
$("#result_main").append("<span style='color:#932192'>" + "> " + "</span>");
$("#result_main").append("<span style='color:blue'>library(QCA)</span><br>");
$("#result_main").append("<span style='color:red'>" + 
    "To cite this package in publications, please use:<br><br>" + 
    "  Dușa, Adrian (2019) QCA with R. A Comprehensive Resource.<br>".split(" ").join("&nbsp;") + 
    "  Springer International Publishing.<br><br>".split(" ").join("&nbsp;"))
createCommandPromptInRconsole();
var keys = {
    cmd: false,
    ctrl: false,
    down: false
};
$("body").on("keydown", function(evt) {
    evt = evt || event; 
    var key = evt.which || evt.keyCode;
    if (key == 17) {
        keys["ctrl"] = true;
    }
    if (key == 91) {
        keys["cmd"] = true;
    }
    if (key == 40) {
        keys["down"] = true;
    }
    if ((keys["ctrl"] || keys["cmd"]) && keys["down"]) {
        if (!spinnerstarted) {
            if(!$("#txtarea").is(":focus")) {
                $("#tempdiv").click();
            }
            $("#result_main").animate({
                scrollTop: $("#result_main")[0].scrollHeight
            }, 1);
            keys = {
                cmd: false,
                ctrl: false,
                down: false
            }
        }
    }
});
$("body").on("keyup", function(evt) {
    evt = evt || event; 
    var key = evt.which || evt.keyCode;
    if (key == 17) {
        keys["ctrl"] = false;
    }
    if (key == 91) {
        keys["cmd"] = false;
    }
    if (key == 40) {
        keys["down"] = false;
    }
});
$("#tempdiv").click();
function showDialogToFront(settings) {
    $("#" + settings.name).show();
    focus = settings.name;
    $("#" + settings.name).css("z-index", "9000");
    var scrollTop = {};
    var scrollLeft = {};
    scrollTop[settings.name + "_main"] = $("#" + settings.name + "_main").scrollTop();
    scrollLeft[settings.name + "_main"] = $("#" + settings.name + "_main").scrollLeft();
    if (settings.inside !== undefined) {
        var keys = getKeys(settings.inside);
        for (var i = 0; i < keys.length; i++) {
            scrollTop[keys[i]] = $("#" + settings.name + "_" + keys[i]).scrollTop();
            scrollLeft[keys[i]] = $("#" + settings.name + "_" + keys[i]).scrollLeft();
        }
    }
    $("#" + settings.name).appendTo(document.body);
    current_command = settings.name;
    if (settings.name != "command" && settings.name != "result") {
        console_command(current_command);
    }
    $("#" + settings.name + "_main").scrollTop(scrollTop[settings.name + "_main"]);
    $("#" + settings.name + "_main").scrollLeft(scrollLeft[settings.name + "_main"]);
    if (settings.inside !== undefined) {
        for (var i = 0; i < keys.length; i++) {
            $("#" + settings.name + "_" + keys[i]).scrollTop(scrollTop[keys[i]]);
            $("#" + settings.name + "_" + keys[i]).scrollLeft(scrollLeft[keys[i]]);
        }
    }
}
function createDialog(dialogsettings) {
    var dialog = document.createElement("div");
    dialog.id = dialogsettings.name;
    document.body.appendChild(dialog);
    current_command = dialogsettings.name;
    $("#" + dialogsettings.name).css({
        position: "absolute",
        width:  dialogsettings.width,
        height: dialogsettings.height
    }).position(dialogsettings.position);
    $("#" + dialogsettings.name).addClass("outerborder");
    if (dialogsettings.resizable) {
        $("#" + dialogsettings.name).resizable({
            resize: function() {
                var twidth = $(this).width();
                var theight = $(this).height();
                if (twidth < dialogsettings.width) {
                    twidth = dialogsettings.width;
                    $("#" + dialogsettings.name).width(twidth);
                }
                if (theight < dialogsettings.height) {
                    theight = dialogsettings.height;
                    $("#" + dialogsettings.name).height(theight);
                }
                showDialogToFront(dialogsettings);
                $("#" + dialogsettings.name + "_main").css({
                    width: twidth + "px",
                    height: (theight - 20) + "px"
                })
                if (dialogsettings.name == "result") {
                    positionSpinner();
                }
            }
        })
    }
    var header = document.createElement("div");
    header.id = dialogsettings.name + "_header";
    dialog = document.getElementById(dialogsettings.name);
    dialog.appendChild(header);
    $("#" + dialogsettings.name + "_header").css({
        fontWeight: "bold",
        padding: "1px 0px 0px 3px"
    }).addClass("header");
    header.innerHTML = dialogsettings.title;
    $("#" + dialogsettings.name + "_header").prepend('<img id="' + dialogsettings.name + '_img' + '" src="css/images/close.png" width="27px"/>')
    $("#" + dialogsettings.name + "_img").css({
        "vertical-align": "middle"
    }).mouseenter(function() {
        $(this).attr('src', 'css/images/closex.png');
    }).mouseleave(function() {
        $(this).attr('src', 'css/images/close.png');
    });
    var clickfired = false;
    $("#" + dialogsettings.name + "_img").mousedown(function(event) {
        clickfired = true;
    });
    $("#" + dialogsettings.name + "_img").click(function(event) {
        event.stopPropagation();
        if (dialogsettings.name == "plotdiv") {
            closeplot += 1;
            plotopen = false;
            Shiny.onInputChange("closeplot", closeplot);
            $("#saveRplot_preview").css({
                "background-image": ""
            });
        }
        if (dialogsettings.name !== "command" && dialogsettings.name !== "result") {
            if (dialogsettings.name == "import") {
                dirfile.filepath = "";
                commobj.read_table.objname = "";
            }
            $("#" + dialogsettings.name).remove();
            if (getKeys(papers).indexOf(dialogsettings.name) >= 0) {
                var dialog_papers = getKeys(papers[dialogsettings.name]);
                for (var i = 0; i < dialog_papers.length; i++) {
                    papers[dialogsettings.name][dialog_papers[i]].remove();
                }
                delete papers[dialogsettings.name];
            }
        }
    });
    $("#" + dialogsettings.name + "_img").mouseup(function(event) {
        if (clickfired) {
            if (dialogsettings.name == "plotdiv") {
                closeplot += 1;
                plotopen = false;
                Shiny.onInputChange("closeplot", closeplot);
                $("#saveRplot_preview").css({
                    "background-image": ""
                });
            }
            if (dialogsettings.name !== "command" && dialogsettings.name !== "result") {
                if (dialogsettings.name == "import") {
                    dirfile.filepath = "";
                    commobj.read_table.objname = "";
                }
                $("#" + dialogsettings.name).remove();
                if (getKeys(papers).indexOf(dialogsettings.name) >= 0) {
                    var dialog_papers = getKeys(papers[dialogsettings.name]);
                    for (var i = 0; i < dialog_papers.length; i++) {
                        papers[dialogsettings.name][dialog_papers[i]].remove();
                    }
                    delete papers[dialogsettings.name];
                }
            }
            clickfired = false;
        }
    });
    var main = document.createElement("div");
    main.id = dialogsettings.name + "_main";
    dialog.appendChild(main);
    $("#" + dialogsettings.name + "_main").css({
        width: ($("#" + dialogsettings.name).width()) + "px",
        height: ($("#" + dialogsettings.name).height() - 20) + "px",
        background: "#ffffff"
    });
    $("#" + dialogsettings.name).draggable({
        handle: "#" + dialogsettings.name + "_header",
        start: function() {
            showDialogToFront(dialogsettings);
        },
        stop: function() {
            if ($(this).offset().top < 31) {
                var left = $(this).offset().left;
                $(this).position({my: "left top", at: "left+" + left + "px top+31px", of: window, collision: "none"});
            }
            var position = $(this).position();
            settings[dialogsettings.name].position.at = "left+" + position.left + "px top+" + position.top + "px";
        }
    });
    if (dialogsettings.resizable && ["data_editor"].indexOf(dialogsettings.name) < 0) {
        $("#" + dialogsettings.name).resizable({
            resize: function() {
                var temph = $(this).height();
                var tempw = $(this).width();
                if (dialogsettings.hasOwnProperty("minHeight")) {
                    if (temph < dialogsettings.minHeight) {
                        $(this).height(dialogsettings.minHeight);
                    }
                }
                if (dialogsettings.hasOwnProperty("minWidth")) {
                    if (tempw < dialogsettings.minWidth) {
                        $(this).width(dialogsettings.minWidth);
                    }
                }
                if (shiftkey) {
                    if (tempw < temph) {
                        $(this).width($(this).height());
                    }
                    if (temph < tempw) {
                        $(this).height($(this).width());
                    }
                }
                if (["command", "result"].indexOf(dialogsettings.name) < 0) {
                    $("#" + dialogsettings.name + "_main").width($(this).width());
                    $("#" + dialogsettings.name + "_main").height($(this).height() - 20);
                }
                else {
                    $("#" + dialogsettings.name + "_main").width($(this).width() - scrollbarsWH);
                    $("#" + dialogsettings.name + "_main").height($(this).height() - 20 - scrollbarsWH);
                }
            }
        })
    }
    if (dialogsettings.inside !== void 0) {
        var keys = getKeys(dialogsettings.inside);
        for (var i = 0; i < keys.length; i++) {
            addDiv(dialogsettings.name, keys[i], dialogsettings.inside[keys[i]]);
        }
    }
    var clientX, clientY;
    $("#" + dialogsettings.name).mousedown(function(event) {
        clientX = event.clientX;
        clientY = event.clientY;
    });
    $("#" + dialogsettings.name).click(function(event) {
        focus = dialogsettings.name;
        if (clientX == event.clientX && clientY == event.clientY) {
            showDialogToFront(dialogsettings);
            if (dialogsettings.name == "result") {
                $("#tempdiv").click();
            }
        }
    });
    $("#" + dialogsettings.name).css("z-index", "9000");
    if (["command", "result"].indexOf(dialogsettings.name) < 0) {
        $("#" + dialogsettings.name).disableTextSelection();
    }
}
$("#result_main").bind("paste", function() {
    if(!$("#txtarea").is(":focus")) {
        var foo = $(this);
        txtcommand += foo.val();
        $("#tempdiv").click();
    }
});
var venn = {
    "s1": [[[500,250,362.5,250,250,362.5,250,500,250,637.5,362.5,750,500,750,637.5,750,750,637.5,750,500,750,362.5,637.5,250,500,250]],[[0],[0]]],
    "s2": [[[500,716,333.5,618.5,333.5,381.5,500,284],[500,716,463.5,738,419,749.5,375,750,237.5,750,125,637.5,125,500,125,362.5,237.5,250,375,250,419,250,463.5,262,500,284],[500,284,536.5,262,581,250,625,250,762.5,250,875,362.5,875,500,875,637.5,762.5,750,625,750,581,749.5,536.5,738,500,716],[500,716,666.5,619,666.5,381.5,500,284]],[[1,2],[2,3],[0,1],[0,3]]],
    "s3": [[[250,391.5,250,256,361,142,500,142,639,142,750,256,750,391.5],[500,824.5,460,847.5,417.5,858,375,858,236.9,858,124.9,746.2,124.9,608.3,124.9,567.2,135.2,527.4,153.1,493,177.4,446.7,212.2,413.8,250,391.5],[250,391.5,250,488.5,304,568,375,608.5],[500,391.5,421.2,346.4,326.8,347.2,250,391.5],[375,608.5,374.5,521,422.5,436.5,500,391.5],[625,608.5,549.5,653,450.4,653,375,608.5],[625,608.5,625.4,695.8,577.3,780.4,500,824.5],[375,608.5,374.6,695.8,422.7,780.4,500,824.5],[625,608.5,625.5,521,577.5,436.5,500,391.5],[500,391.5,578.8,346.4,673.2,347.2,750,391.5],[750,391.5,750,488.5,696,568,625,608.5],[750,391.5,787.8,413.8,822.6,446.7,846.9,493,864.8,527.4,875.1,567.2,875.1,608.3,875.1,746.2,763.1,858,625,858,582.5,858,540,847.5,500,824.5]],[[0,11,1],[6,11,10],[0,9,3],[8,9,10],[1,2,7],[5,7,6],[2,4,3],[4,8,5]]],
    "s4": [[[359.9,344.3,406.8,297.3,453.7,250.4,500.8,203.3],[641.7,344.2,594.7,297.3,547.9,250.3,500.8,203.3],[765.2,467.7,724.1,426.6,682.9,385.5,641.7,344.2],[236.5,467.7,283,514.3,329.7,560.9,376.4,607.6],[625.3,607.7,583.8,566.2,542.3,524.7,500.8,483.2],[376.4,607.6,417.9,566.1,459.4,524.7,500.8,483.2],[359.9,344.3,365.8,348.9,371.5,353.9,376.9,359.3,418.7,400.8,459.5,441.9,500.8,483.2],[500.8,483.2,542,441.9,583.3,400.7,624.5,359.4,630,353.9,635.7,348.8,641.7,344.2],[691.2,769.8,711.7,761.2,730.9,748.5,747.7,731.7,753.4,726,759.1,720.4,764.8,714.7,833.3,646.1,833.7,536.2,765.2,467.7],[691.2,769.8,697.9,717.9,681.3,663.6,641.2,623.6,635.9,618.3,630.6,613,625.3,607.7],[310.3,770.2,303.5,718.2,320.2,663.8,360.3,623.7,365.7,618.3,371,613,376.4,607.6],[500.8,939.5,456.1,939.5,411.5,922.1,377.1,887.7,371.5,882.2,366.2,876.7,360.5,871.1,331.9,842.9,315.1,807.1,310.3,770.2],[500.8,939.5,545.5,939.6,590.7,922.8,624.6,888,630.1,882.3,635.9,876.9,641.4,871.2,670,842.7,686.4,806.7,691.2,769.8],[310.3,770.2,273.5,753.6,263.6,742.1,236.5,715,168,646.5,168,536.2,236.5,467.7],[186.3,321.2,194.2,270.3,213.4,240.3,253.1,203,321.6,134.5,432.3,134.8,500.8,203.3],[500.8,203.3,569.3,134.8,679.6,134.8,748.2,203.3,791.3,244.1,806.9,272.5,815.1,321.4],[815.1,321.4,845,335.1,854.7,341.7,888.9,376.5,957.4,445,957.4,555.3,888.9,623.8,800.8,711.9,712.7,799.9,624.6,888,590.3,922.4,545.5,939.6,500.8,939.5],[500.8,939.5,456.1,939.5,411.5,922.1,377.1,887.7,289,799.6,201.9,712.5,112.8,623.4,43.17,554.9,43.17,444.6,112.2,376.7,134.3,354.5,149.9,336.9,186.3,321.2],[765.2,467.7,718.6,514.3,672,561,625.3,607.7],[376.4,607.6,418,649.2,459.4,690.7,500.8,732.1],[625.3,607.7,583.8,649.2,542.3,690.6,500.8,732.1],[236.5,467.7,277.6,426.6,318.6,385.5,359.9,344.3],[641.7,344.2,692.1,305.4,758.4,297.8,815.1,321.4],[359.9,344.3,309.6,305.4,242.9,297.8,186.3,321.2],[815.1,321.4,821.9,373.3,805.2,427.6,765.2,467.7],[500.8,732.1,449.1,783.9,373.7,796.5,310.3,770.2],[186.3,321.2,179.7,373.2,196.5,427.7,236.5,467.7],[500.8,732.1,552.5,783.8,627.9,796.3,691.2,769.8]],[[14,15,16,17],[8,24,16,12],[1,15,22],[2,22,24],[0,14,23],[8,18,9],[0,1,7,6],[2,7,4,18],[11,13,26,17],[11,25,27,12],[3,10,13],[10,19,25],[21,23,26],[9,20,27],[3,5,6,21],[4,5,19,20]]],
    "s5": [[[757,276.8,705.7,268.3,668.3,268.6,620.7,286.6],[620.7,286.6,590.6,297.2,567.3,309.2,548.2,322.1],[548.2,322.1,516.4,343.7,496.2,367.8,475.3,392.1],[475.3,392.1,443.3,431.9,407.9,457.6,367.2,481.1],[367.2,481.1,342.6,495.2,316.9,509.4,295.2,525.2],[295.2,525.2,265,547.3,242.7,572.6,242.3,606.1],[242.3,606.1,241.9,654.1,274.2,686.8,325,693.5],[325,693.5,360.3,698.1,402.8,693,443.9,690.5],[443.9,690.5,477,688.5,507.7,689.4,538,695.7],[538,695.7,566.2,701.4,590.1,709.8,615.7,716],[615.7,716,693.3,730,732.5,723.8,787.8,701.9],[757,276.8,831.6,294.8,894.8,350.7,918.6,429.2,952.5,540.9,894.6,659.2,787.8,701.9],[348.9,169.4,321,211.9,309.2,247.4,310.4,307.6],[310.4,307.6,311.2,340.9,315.5,367.7,322.2,390.5],[322.2,390.5,333.1,427.7,350,454.3,367.2,481.1],[367.2,481.1,395.6,522.6,411.6,560,422.5,605.3],[422.5,605.3,429.3,633.7,435.2,663.7,443.9,690.5],[443.9,690.5,455.1,725.4,470.9,754.7,499.5,767.5],[499.5,767.5,553.4,788,592.4,760.5,615.7,716],[615.7,716,630,683.9,637.4,642.9,646.6,603.8],[646.6,603.8,654.4,569.3,665.5,538.2,681.4,509.1],[681.4,509.1,693.5,487.5,706.3,469.3,718,449.8],[718,449.8,748.1,395.6,761.6,364.3,757,276.8],[348.9,169.4,389.2,110.3,457.1,71.49,533.9,71.49,651.2,71.49,747.8,161.9,757,276.8],[393.5,855.2,440.8,832.4,472.6,805.3,499.5,767.5],[499.5,767.5,518.2,741.1,530.1,717.7,538,695.7],[538,695.7,551,659.9,553.2,628.3,556.4,596],[556.4,596,559.5,546,574.1,504.5,594,462.8],[594,462.8,606,437.6,618.9,412,627.5,387.2],[627.5,387.2,640.6,349.9,644,314.9,620.7,286.6],[620.7,286.6,580.5,238.8,540.9,249.5,505.2,268.6],[505.2,268.6,474.2,284.9,443.1,312.2,411.4,337],[411.4,337,383.1,359.6,353.7,377.7,322.2,390.5],[322.2,390.5,298.9,399.9,277.2,406.3,255.3,415],[255.3,415,191.6,450.3,156.6,466.8,120.5,527.6],[393.5,855.2,321.3,884.4,235.7,875,170.1,823.4,79.32,750.8,59.69,622.6,120.5,527.6],[120.5,527.6,158.3,568,195.7,592.2,242.3,606.1],[242.3,606.1,272.2,615.1,297.5,618.9,320.3,619.7],[320.3,619.7,359.4,621,390.4,612.9,422.5,605.3],[422.5,605.3,469.7,592.3,511.4,591.1,556.4,596],[556.4,596,586.6,599.3,617.9,603.9,646.6,603.8],[646.6,603.8,684.3,603.7,717.3,595.8,736.8,566.7],[736.8,566.7,766.2,522.8,756.5,482.3,718,449.8],[718,449.8,693.6,427.3,659.7,407.7,627.5,387.2],[627.5,387.2,596.6,367.5,571,347.8,548.2,322.1],[548.2,322.1,531.9,303.9,519.2,285.9,505.2,268.6],[505.2,268.6,461.3,214.8,422.8,187.6,348.9,169.4],[120.5,527.6,74.02,469,57.84,388.6,84.78,313.3,123.4,203.4,238.9,143.2,348.9,169.4],[787.8,701.9,781.1,647.5,765.9,604.4,736.8,566.7],[736.8,566.7,717.8,541.5,699.3,523.1,681.4,509.1],[681.4,509.1,651.7,485.9,623,474.3,594,462.8],[594,462.8,545.9,445.3,510.7,423.4,475.3,392.1],[475.3,392.1,454.2,373.5,433,353.3,411.4,337],[411.4,337,378.5,312,345,295.9,310.4,307.6],[310.4,307.6,263,327.7,243.7,363.3,255.3,415],[255.3,415,262.4,449.5,280.2,487.6,295.2,525.2],[295.2,525.2,308.1,557,316.8,587.8,320.3,619.7],[320.3,619.7,323.1,645.8,323,669.4,325,693.5],[325,693.5,340.2,765.1,348.1,799.5,393.5,855.2],[787.8,701.9,790.5,775.2,757.2,848.2,693,893.1,597.7,960.7,467.4,942.6,393.5,855.2]],[[11,59,35,47,23],[6,58,35,36],[10,59,24,18],[7,17,24,58],[11,48,42,22],[6,57,37],[10,48,41,19],[7,16,38,57],[0,30,46,23],[1,45,30],[9,18,25],[8,25,17],[0,29,43,22],[1,44,29],[9,19,40,26],[8,26,39,16],[12,54,34,47],[5,36,34,55],[13,33,54],[4,55,33,14],[21,42,49],[5,37,56],[20,49,41],[4,56,38,15],[12,53,31,46],[2,52,31,45],[13,32,53],[3,14,32,52],[21,43,28,50],[2,51,28,44],[20,50,27,40],[3,15,39,27,51]]],
    "s6": [[[556,151,610,127,670,135,685,192],[685,192,689,206,692,222,693,237],[693,237,695,263,695,276,695,302],[695,302,695,333,695,364,695,406],[695,406,694,462,685,592,664,680],[664,680,656,713,646,734,631,757],[631,757,614,782,600,799,580,817],[580,817,558,837,528,855,497,869],[497,869,457,889,399,897,361,898],[361,898,196,898,51,766,81,561],[81,561,87,522,98,493,116,466],[116,466,129,447,145,433,162,420],[162,420,176,409,195,397,210,387],[210,387,254,359,270,337,275,297],[275,297,278,277,271,252,285,222],[285,222,305,182,342,169,385,182],[385,182,432,196,459,191,495,177],[495,177,516,168,535,160,556,151],[767,211,918,266,1022,444,922,611],[922,611,864,700,826,735,748,759],[748,759,719,765,703,766,673,764],[673,764,654,762,644,760,631,757],[631,757,568,742,555,759,523,795],[523,795,506,815,481,839,464,845],[464,845,428,860,391,853,366,818],[366,818,348,791,337,773,291,768],[291,768,262,765,231,762,202,754],[202,754,146,738,119,701,136,642],[136,642,143,619,149,609,158,592],[158,592,173,567,176,560,190,540],[190,540,204,520,224,493,232,482],[232,482,306,389,391,305,455,255],[455,255,478,237,488,230,522,214],[522,214,542,205,563,198,585,194],[585,194,623,189,653,188,685,192],[685,192,718,196,745,203,767,211],[747,312,765,316,785,321,803,327],[803,327,857,347,867,398,849,436],[849,436,840,454,830,472,813,492],[813,492,796,511,788,520,767,543],[767,543,739,577,726,595,733,646],[733,646,736,667,740,676,743,698],[743,698,746,719,748,741,748,759],[748,759,748,796,740,822,723,849],[723,849,657,966,453,977,361,898],[361,898,337,879,323,857,313,838],[313,838,302,818,295,801,291,768],[291,768,289,748,288,742,288,722],[288,722,288,672,270,667,226,639],[226,639,207,627,177,609,158,592],[158,592,136,573,127,543,151,506],[151,506,173,472,171,448,162,420],[162,420,152,387,146,361,145,327],[145,327,143,265,168,209,240,215],[240,215,249,216,269,219,285,222],[285,222,306,226,336,232,352,235],[352,235,382,241,421,248,455,255],[455,255,511,266,552,274,625,288],[625,288,645,292,670,297,695,302],[695,302,713,306,728,308,747,312],[695,406,705,414,746,444,763,455],[763,455,774,463,791,475,813,492],[813,492,828,504,847,519,864,536],[864,536,885,557,908,583,922,611],[922,611,981,733,865,862,723,849],[723,849,692,847,675,844,645,837],[645,837,620,831,595,823,580,817],[580,817,557,808,541,802,523,795],[523,795,480,780,458,776,412,796],[412,796,401,801,386,808,366,818],[366,818,357,823,334,833,313,838],[313,838,241,853,202,803,202,754],[202,754,202,736,206,715,211,695],[211,695,216,676,221,658,226,639],[226,639,236,599,228,570,190,540],[190,540,180,531,163,517,151,506],[151,506,136,493,125,480,116,466],[116,466,86,416,99,354,145,327],[145,327,162,318,183,313,209,311],[209,311,237,309,257,304,275,297],[275,297,318,281,336,260,352,235],[352,235,363,218,374,198,385,182],[385,182,399,160,412,143,435,130],[435,130,474,110,526,117,556,151],[556,151,568,164,575,176,585,194],[585,194,593,208,601,230,608,245],[608,245,614,258,619,275,625,288],[625,288,649,343,654,373,695,406],[608,245,632,243,669,238,693,237],[693,237,737,236,746,259,747,312],[747,312,748,342,754,373,792,398],[792,398,809,410,831,424,849,436],[849,436,897,467,894,497,864,536],[864,536,846,559,824,580,800,600],[800,600,787,612,757,633,733,646],[733,646,711,658,689,669,664,680],[664,680,571,723,456,742,361,734],[361,734,328,731,315,728,288,722],[288,722,269,717,232,705,211,695],[211,695,179,681,158,663,136,642],[136,642,105,608,97,593,81,561],[81,561,-5,388,68,126,320,117],[320,117,380,117,411,122,435,130],[435,130,456,138,479,155,495,177],[495,177,503,188,513,203,522,214],[522,214,540,236,562,250,608,245],[803,327,805,354,803,372,792,398],[792,398,783,419,775,433,763,455],[763,455,746,487,745,503,767,543],[767,543,778,564,790,577,800,600],[800,600,820,642,800,687,743,698],[743,698,715,704,685,710,673,764],[673,764,668,788,658,815,645,837],[645,837,626,874,565,908,497,869],[497,869,483,860,473,852,464,845],[464,845,448,833,428,812,412,796],[412,796,394,777,377,757,361,734],[361,734,306,655,265,577,232,482],[232,482,223,456,214,415,210,387],[210,387,207,365,207,335,209,311],[209,311,212,277,223,244,240,215],[240,215,260,176,288,142,320,117],[320,117,448,16,683,41,767,211],[767,211,781,239,800,282,803,327]],[[9,101,122,18,64,44],[9,100,27,71,45],[8,44,65,113],[8,45,70,24,114],[19,43,64],[26,71,46],[20,112,65,43],[25,46,70],[18,63,92,37,123],[27,99,72],[37,91,106],[24,69,115],[19,42,110,93,63],[26,72,98,47],[20,111,42],[25,47,97,116,69],[0,35,122,102,83],[0,34,84],[7,113,66],[7,114,23,67],[17,83,103],[17,84,33,104],[6,66,112,21],[6,67,22],[1,89,36,123,35],[1,88,85,34],[36,106,90],[23,115,68],[41,110,94],[33,85,105],[5,21,111,41,95],[5,22,68,116,96],[10,77,53,121,101],[10,76,50,28,100],[53,120,78],[29,75,50],[11,52,77],[11,51,76],[12,119,78,52],[12,118,30,75,51],[38,62,92],[28,49,73,99],[38,61,107,91],[29,74,49],[39,109,93,62],[48,73,98],[39,108,61],[30,117,97,48,74],[15,82,102,121,54],[15,81,55],[14,54,120,79],[14,55,80],[16,103,82],[16,104,32,56,81],[13,79,119],[13,80,56,31,118],[2,59,89],[2,58,86,88],[3,60,107,90,59],[3,87,58],[40,94,109],[32,105,86,57],[4,95,40,108,60],[4,96,117,31,57,87]]],
    "s7": [[[179,584,197,592,213,595,230,600],[230,600,245,604,258,606,272,607],[272,607,288,608,297,608,316,609],[316,609,351,610,370,630,368,668],[368,668,367,685,366,698,366,713],[366,713,366,730,366,745,367,760],[367,760,368,776,371,793,375,808],[375,808,379,826,383,837,389,853],[718,812,684,1014,450,1016,389,853],[718,812,721,795,721,776,721,760],[721,760,721,742,720,726,717,709],[717,709,714,693,713,682,709,666],[709,666,706,644,712,633,736,628],[736,628,752,625,765,625,781,623],[781,623,815,618,819,595,797,571],[797,571,788,560,776,551,768,540],[768,540,756,521,757,508,772,493],[772,493,783,483,796,476,807,467],[807,467,837,443,826,424,804,413],[804,413,790,406,778,402,760,399],[760,399,740,397,726,397,711,397],[711,397,694,398,681,398,668,399],[668,399,603,402,556,377,534,312],[534,312,528,295,524,284,517,269],[517,269,510,255,504,243,496,230],[496,230,488,217,481,205,472,194],[472,194,460,180,449,167,432,155],[432,155,387,123,331,140,322,205],[322,205,320,226,320,240,321,262],[321,262,322,278,324,291,327,306],[327,306,331,321,334,332,338,346],[338,346,349,384,346,405,304,415],[304,415,288,418,272,421,257,425],[257,425,241,429,222,435,207,441],[207,441,193,447,183,451,168,461],[168,461,116,494,118,561,179,584],[496,153,486,167,483,171,472,194],[472,194,467,205,461,216,453,238],[453,238,448,253,446,262,441,278],[441,278,429,320,417,328,372,317],[372,317,355,313,342,310,327,306],[327,306,311,302,296,299,280,297],[280,297,265,295,250,294,236,294],[236,294,219,295,203,296,185,299],[150,627,-16,568,-9,333,185,299],[150,627,168,634,185,638,204,642],[204,642,220,646,235,648,251,649],[251,649,266,650,280,650,294,651],[294,651,317,652,327,660,327,683],[327,683,327,700,324,713,323,729],[323,729,321,761,341,774,367,760],[367,760,382,751,396,743,408,736],[408,736,427,725,441,728,450,748],[450,748,455,759,461,776,468,789],[468,789,482,816,503,821,525,794],[525,794,533,784,540,775,547,756],[547,756,554,739,557,728,560,714],[560,714,563,699,565,681,568,667],[568,667,581,604,624,571,685,560],[685,560,701,557,715,554,730,551],[730,551,745,548,754,545,768,540],[768,540,785,533,798,527,812,520],[812,520,828,512,838,505,851,495],[851,495,905,453,893,395,832,376],[832,376,814,370,803,367,780,364],[780,364,754,361,747,361,730,361],[730,361,716,361,703,362,689,362],[689,362,645,364,629,345,633,307],[633,307,635,289,635,275,635,259],[635,259,635,241,634,226,631,212],[631,212,626,188,621,177,616,165],[616,165,594,108,535,100,496,153],[763,271,748,273,738,273,718,278],[718,278,703,281,689,286,675,291],[675,291,661,296,646,301,633,307],[633,307,600,321,583,316,563,274],[563,274,552,249,551,245,545,232],[545,232,537,216,528,200,522,190],[522,190,508,169,504,163,496,153],[496,153,475,129,468,121,456,110],[185,299,121,101,331,-10,456,110],[185,299,190,314,196,331,203,344],[203,344,212,361,219,372,231,390],[231,390,239,402,248,413,257,425],[257,425,272,446,268,458,243,472],[243,472,231,479,220,485,209,493],[209,493,185,510,185,530,213,546],[213,546,230,555,241,558,258,565],[258,565,282,575,282,586,272,607],[272,607,264,622,257,635,251,649],[251,649,241,675,251,690,281,690],[281,690,296,690,315,686,327,683],[327,683,341,679,355,674,368,668],[368,668,378,662,396,655,408,648],[408,648,465,620,525,625,568,667],[568,667,580,679,589,687,600,697],[600,697,611,706,620,714,632,723],[632,723,647,733,664,742,678,748],[678,748,692,753,706,757,721,760],[721,760,782,771,823,731,802,667],[802,667,796,648,788,635,781,623],[781,623,770,605,760,592,753,582],[753,582,746,571,739,562,730,551],[730,551,705,521,705,494,734,473],[734,473,747,464,757,458,770,447],[770,447,785,434,793,426,804,413],[804,413,816,400,822,391,832,376],[832,376,870,320,842,265,763,271],[845,550,828,532,824,529,812,520],[812,520,796,507,785,500,772,493],[772,493,757,486,748,481,734,473],[734,473,699,454,691,434,711,397],[711,397,722,378,724,373,730,361],[730,361,737,347,744,332,749,318],[749,318,755,302,759,289,763,271],[763,271,767,250,769,237,771,214],[456,110,545,-35,791,29,771,214],[456,110,448,123,441,133,432,155],[432,155,426,170,420,187,415,206],[415,206,411,222,409,233,406,249],[406,249,400,278,394,285,364,276],[364,276,346,271,335,266,321,262],[321,262,287,251,272,264,280,297],[280,297,285,313,287,317,296,339],[296,339,307,365,300,373,277,380],[277,380,258,385,247,387,231,390],[231,390,196,396,188,414,207,441],[207,441,214,451,228,462,243,472],[243,472,256,481,268,486,280,492],[280,492,295,499,310,507,322,513],[322,513,375,541,412,572,408,648],[408,648,407,662,406,673,405,692],[405,692,405,705,404,712,408,736],[408,736,410,751,413,766,417,780],[417,780,422,796,426,807,433,822],[433,822,459,876,529,888,561,831],[561,831,570,816,575,802,580,788],[580,788,587,768,588,764,592,744],[592,744,596,723,599,709,600,697],[600,697,605,656,633,646,667,655],[667,655,682,658,693,662,709,666],[709,666,719,669,740,671,756,671],[756,671,774,671,785,670,802,667],[802,667,877,653,889,600,845,550],[673,796,676,778,677,766,678,748],[678,748,679,731,678,714,675,697],[675,697,673,682,671,669,667,655],[667,655,658,619,669,598,708,589],[708,589,724,586,739,584,753,582],[753,582,768,580,782,576,797,571],[797,571,813,566,831,557,845,550],[845,550,859,542,872,535,885,526],[771,214,965,190,1050,415,885,526],[771,214,754,216,739,221,723,225],[723,225,706,229,691,234,676,240],[676,240,661,247,649,253,635,259],[635,259,610,270,602,265,592,242],[592,242,586,225,581,211,576,199],[576,199,566,173,541,166,522,190],[522,190,511,205,508,211,496,230],[496,230,483,249,472,251,453,238],[453,238,437,225,424,215,415,206],[415,206,391,182,363,192,362,221],[362,221,361,238,361,241,364,276],[364,276,367,295,369,303,372,317],[372,317,375,329,376,334,381,358],[381,358,401,440,384,470,322,513],[322,513,311,521,301,527,291,535],[291,535,283,541,271,552,258,565],[258,565,243,580,238,589,230,600],[230,600,222,610,212,626,204,642],[204,642,177,694,209,747,272,740],[272,740,289,738,304,735,323,729],[323,729,339,724,352,719,366,713],[366,713,380,706,392,699,405,692],[405,692,434,673,469,677,483,709],[483,709,491,727,494,736,500,752],[500,752,503,760,513,781,525,794],[525,794,540,813,545,818,561,831],[561,831,609,870,662,861,673,796],[375,808,399,796,405,790,417,780],[417,780,430,770,445,754,450,748],[450,748,461,734,472,721,483,709],[483,709,506,686,537,689,560,714],[560,714,571,726,581,735,592,744],[592,744,602,753,616,763,626,770],[626,770,641,781,655,788,673,796],[673,796,687,802,704,808,718,812],[885,526,1025,670,895,870,718,812],[885,526,873,515,864,505,851,495],[851,495,836,484,822,476,807,467],[807,467,794,460,782,454,770,447],[770,447,752,434,749,421,760,399],[760,399,767,386,774,375,780,364],[780,364,797,335,780,316,749,318],[749,318,735,319,721,321,707,321],[707,321,685,322,675,314,675,291],[675,291,676,269,677,258,676,240],[676,240,674,208,661,197,631,212],[631,212,617,220,603,231,592,242],[592,242,576,259,573,263,563,274],[563,274,553,287,543,300,534,312],[534,312,495,365,450,379,381,358],[381,358,367,353,351,349,338,346],[338,346,325,343,314,341,296,339],[296,339,284,338,266,337,253,338],[253,338,235,339,218,341,203,344],[203,344,135,359,121,420,168,461],[168,461,181,473,189,480,209,493],[209,493,217,498,233,508,247,515],[247,515,257,520,273,527,291,535],[291,535,324,550,336,577,316,609],[316,609,307,624,301,634,294,651],[294,651,288,666,284,680,281,690],[281,690,276,707,273,723,272,740],[272,740,268,806,316,839,375,808],[236,294,241,309,246,321,253,338],[253,338,257,348,264,362,277,380],[277,380,287,393,296,403,304,415],[304,415,328,451,320,466,280,492],[280,492,269,499,258,506,247,515],[247,515,235,525,224,536,213,546],[213,546,198,561,188,572,179,584],[179,584,169,598,159,611,150,627],[389,853,214,959,50,785,150,627],[389,853,405,843,418,834,433,822],[433,822,445,811,457,801,468,789],[468,789,479,776,488,765,500,752],[500,752,515,736,530,736,547,756],[547,756,559,769,568,780,580,788],[580,788,600,802,620,797,626,770],[626,770,630,755,629,743,632,723],[632,723,635,700,653,690,675,697],[675,697,688,700,703,706,717,709],[717,709,754,717,763,697,756,671],[756,671,751,652,743,637,736,628],[736,628,723,609,718,602,708,589],[708,589,699,578,694,572,685,560],[685,560,649,509,639,456,668,399],[668,399,675,387,682,374,689,362],[689,362,697,348,702,335,707,321],[707,321,714,302,716,289,718,278],[718,278,722,258,722,244,723,225],[723,225,726,161,668,129,616,165],[616,165,605,173,592,183,576,199],[576,199,565,210,554,221,545,232],[545,232,534,245,527,256,517,269],[517,269,494,301,470,302,441,278],[441,278,419,258,420,260,406,249],[406,249,394,240,375,227,362,221],[362,221,345,213,337,210,322,205],[322,205,253,184,212,221,236,294]],[[8,188,152,116,80,44,224],[35,223,44,81,207],[7,224,45,171,215],[0,170,45,223],[8,187,179,135,225],[35,222,86,208],[7,225,134,180],[0,169,87,222],[9,99,143,151,188],[34,207,82,126],[6,215,172,50],[5,50,173],[9,98,144,187],[34,208,85,127],[6,180,133,51],[5,51,132,174],[63,107,115,152,189],[18,106,63,190],[46,90,214,171],[1,89,46,170],[136,230,186,179],[18,105,191],[137,185,230],[1,88,169],[62,189,151,108],[17,190,62,109],[49,172,214,91],[4,173,49,92],[97,144,186,231],[17,191,104,110],[96,231,185,138],[4,174,131,93],[71,79,116,153,243],[71,78,158,244],[70,243,154,198],[70,244,157,199],[54,178,135,226],[54,177,227],[53,226,134,181],[53,227,176,182],[14,150,143,100],[33,126,83],[69,198,155],[69,199,156],[14,149,101],[33,127,84],[52,181,133],[52,182,175,132],[72,242,153,115],[77,158,245],[73,197,154,242],[76,245,157,200],[55,229,136,178],[55,228,177],[56,184,137,229],[56,183,176,228],[15,61,108,150],[16,109,61],[68,155,197,74],[68,156,200,75],[15,60,102,149],[16,110,103,60],[57,95,138,184],[57,94,131,175,183],[27,251,43,80,117],[43,81,206,216],[28,122,42,251],[42,216,205,123],[27,250,162,118],[86,221,209],[28,121,163,250],[87,168,210,221],[10,234,142,99],[82,125,217,206],[11,141,234],[124,217,205],[10,233,145,98],[85,209,220,128],[11,140,146,233],[129,167,210,220],[64,194,114,107],[19,193,64,106],[47,213,90],[2,212,47,89],[119,249,162],[19,192,105],[120,163,249],[2,211,168,88],[65,113,194],[20,112,65,193],[48,91,213],[3,92,48,212],[97,145,232],[20,111,104,192],[96,232,146,139],[3,93,130,167,211],[26,117,79,36],[25,36,78,159],[29,41,122],[30,204,123,41],[26,118,161,37],[25,37,160],[29,40,164,121],[30,203,165,40],[13,100,142,235],[32,83,125,218],[12,235,141],[31,218,124,204],[13,101,148,236],[32,84,128,219],[12,236,147,140],[31,219,129,166,203],[72,241,195,114],[24,159,77,246],[73,196,241],[23,246,76,201],[38,248,119,161],[24,160,38,247],[39,164,120,248],[23,247,39,165,202],[66,240,195,113],[21,239,66,112],[67,74,196,240],[22,201,75,67,239],[59,102,148,237],[21,238,59,103,111],[58,237,147,139,95],[22,202,166,130,94,58,238]]],
    "l1": {
        "x": [485],
        "y": [213]
    },
    "l2": {
        "x": [100, 888],
        "y": [300, 300]
    },
    "l3": {
        "x": [100, 500, 900],
        "y": [440,  90, 440]
    },
    "l4": {
        "x": [ 75, 220, 780, 925],
        "y": [330, 140, 140, 330]
    },
    "l5": {
        "x": [ 90, 533, 885, 700, 163],
        "y": [200,  38, 280, 960, 890]
    },
    "l6": {
        "x": [100, 500, 910, 925, 550, 100],
        "y": [140,  25, 225, 835, 970, 860]
    },
    "l7": {
        "x": [220, 685, 935, 935, 600, 155,  50],
        "y": [ 45,  23, 220, 800, 985, 880, 310]
    }
}
}); 
