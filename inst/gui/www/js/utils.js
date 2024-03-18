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

navigator.browserType = (function(){
    var N = navigator.appName, ua = navigator.userAgent, tem;
    var M = ua.match(/(opera|chrome|safari|firefox|msie|trident)\/?\s*(\.?\d+(\.\d+)*)/i);
    if (M && (tem= ua.match(/version\/([\.\d]+)/i))!= null) {M[2]=tem[1];}
    M = M? M[1] : N;
    return M;
})();
$.extend($.fn.disableTextSelection = function() {
    return this
         .attr('unselectable', 'on')
         .css('user-select', 'none')
         .on('selectstart', false);
});
function makePapers(obj, papers) {
    papers[obj.name] = new Array();
    papers[obj.name]["main"] = Raphael(obj.name + "_main", obj.width, obj.height);
    if (obj.inside !== undefined) {
        var keys = getKeys(obj.inside);
        for (var i = 0; i < keys.length; i++) {
            papers[obj.name][keys[i]] = Raphael(obj.name + "_" + keys[i], obj.inside[keys[i]].width, obj.inside[keys[i]].height);
        }
    }
}
function addDiv(parent, child, settings) {
    var div1 = document.getElementById(parent + "_main");
    var div2 = document.createElement("div");
    div2.id = parent + "_" + child + ((settings.border)?"_border":"");
    div1.appendChild(div2);
    $("#" + div2.id).css({
        position: "absolute",
        left:     (settings.left   + "px"),
        top:      (settings.top    + "px"),
        width:    (settings.width  + "px"),
        height:   (settings.height + "px")
    });
    if (settings.border) {
        var div3 = document.getElementById(parent + "_" + child + "_border");
        var div4 = document.createElement("div");
        div4.id = parent + "_" + child;
        div3.appendChild(div4);
        $("#" + parent + "_" + child).css({
            "max-height": (settings.height + "px"),
            "overflow-x": "hidden",
            "overflow-y": ((child == "path")?"hidden":"auto")
        });
        if (child != "path" && child != "direxp") {
            $("#" + parent + "_" + child + "_border").addClass("border");
        }
    }
}
Raphael.fn.counter = function(options) {
    if (options.fontsize == void 0) {
        options.fontsize = 14;
    }
    if (options.width == void 0) {
        options.width = 18;
    }
    var cntr = new Array();
    cntr.active = true;
    cntr.value = options.startval;
    var txtanchor = "middle";
    cntr.textlabel = this.text(options.x, options.y, "")
        .attr({"text-anchor": txtanchor, "font-size": options.fontsize + "px"});
    cntr.textvalue = this.text(options.x, options.y, "" + options.startval)
        .attr({"text-anchor": txtanchor, "font-size": options.fontsize + "px"});
    cntr.downsign = this.path([
        ["M", options.x - 12 - options.width / 2, options.y - options.textheight/4],
        ["l", 12, 0],
        ["l", -6, 12],
        ["z"]
    ]).attr({fill: "#eeeeee", "stroke-width": 1.2, stroke: "#a0a0a0"});
    cntr.upsign = this.path([
        ["M", options.x + options.width / 2, options.y + options.textheight/4],
        ["l", 12, 0],
        ["l", -6, -12],
        ["z"]
    ]).attr({fill: "#eeeeee", "stroke-width": 1.2, stroke: "#a0a0a0"});
    cntr.down = this.rect(options.x - 22, options.y - 6, 15, 15)
        .attr({fill: "#fff", opacity: 0, stroke: "#000", "stroke-width": 1, cursor: "pointer"})
        .click(function() {
            if (cntr.value > options.startval) {
                cntr.value -= 1;
                cntr.textvalue.attr({"text": ("" + cntr.value)});
            }
        });
    cntr.up = this.rect(options.x + 8, options.y - 8, 15, 15)
        .attr({fill: "#fff", opacity: 0, stroke: "#000", "stroke-width": 1, cursor: "pointer"})
        .click(function() {
            if (cntr.value < options.maxval) {
                cntr.value += 1;
                cntr.textvalue.attr({"text": ("" + cntr.value)});
            }
        });
    cntr.hideIt = function() {
        cntr.upsign.hide();
        cntr.downsign.hide();
        cntr.up.hide();
        cntr.down.hide();
        cntr.textvalue.hide();
        cntr.textlabel.hide();
    }
    cntr.showIt = function() {
        cntr.upsign.show();
        cntr.downsign.show();
        cntr.up.show();
        cntr.down.show();
        cntr.textvalue.show();
        cntr.textlabel.show();
    }
    cntr.label = function(options) {
        if (options.anchor == void 0) {
            options.anchor = "end";
        }
        if (options.label != void 0) {
            cntr.textlabel.attr({"text": "" + options.label, "text-anchor": options.anchor});
        }
        if (options.x != void 0 && options.y != void 0) {
            cntr.textlabel.transform("t" + options.x + "," + options.y);
        }
    }
    return(cntr);
}
Raphael.fn.checkBox = function(options) {
    if (missing(options.dim)) {
        options.dim = 12;
    }
    if (missing(options.pos)) {
        options.pos = 3;
    }
    if (missing(options.fontsize)) {
        options.fontsize = 14;
    }
    var cb = new Array();
    cb.active = true;
    cb.label = new Array(1);
    var txtanchor = "start";
    var xpos = options.x;
    var ypos = options.y;
    if (options.pos == 1) { 
        xpos -= 8;
        ypos += options.dim / 2;
        txtanchor = "end";
    }
    else if (options.pos == 2) { 
        xpos += options.dim / 2;
        ypos -= options.dim;
        txtanchor = "middle";
    }
    else if (options.pos == 3) { 
        xpos += 20;
        ypos += options.dim / 2;
    }
    else { 
        xpos += options.dim / 2;
        ypos += 5;
        txtanchor = "middle";
    }
    cb.label[0] = this.text(xpos, ypos, options.label)
        .attr({"text-anchor": txtanchor, "font-size": (options.fontsize + "px")});
    cb.box = this.rect(options.x, options.y, options.dim, options.dim)
        .attr({fill: options.isChecked?"#97bd6c":"#eeeeee","stroke-width": 1.2, stroke: "#a0a0a0"});
    cb.chk = this.path([
        ["M", options.x + 0.2*options.dim, options.y + 0.3*options.dim],
        ["l", 0.15*2*options.dim, 0.2*2*options.dim],
        ["l", 0.3*2*options.dim, -0.45*2*options.dim]
    ]).attr({"stroke-width": 2});
    if (options.isChecked) {
        cb.box.attr({fill: "#97bd6c"});
        cb.chk.show();
    }
    else {
        cb.box.attr({fill: "#eeeeee"});
        cb.chk.hide();
    }
    cb.isChecked = options.isChecked;
    cb.cover = this.rect(options.x, options.y, options.dim, options.dim)
        .attr({fill: "#fff", opacity: 0, cursor: "pointer"})
        .click(function() {
            if (cb.active) {
                cb.isChecked = !cb.isChecked;
                this.isChecked = cb.isChecked;
                if (cb.isChecked) {
                    cb.box.attr({fill: "#97bd6c"});
                    cb.chk.show();
                }
                else {
                     cb.box.attr({fill: "#eeeeee"});
                     cb.chk.hide();
                }
            }
        });
    cb.cover.active = true;
    cb.activate = function() {
        cb.active = true;
        cb.cover.active = true;
        cb.cover.attr({fill: "#000", opacity: 0, cursor: "pointer"});
    }
    cb.deactivate = function() {
        cb.active = false;
        cb.cover.active = false;
        cb.cover.attr({fill: "#000", opacity: 0.2, cursor: "default"});
    }
    cb.uncheck = function() {
        cb.isChecked = false;
        cb.box.attr({fill: "#eeeeee"});
        cb.chk.hide();
        cb.cover.isChecked = false;
    }
    cb.check = function() {
        cb.isChecked = true;
        cb.box.attr({fill: "#97bd6c"});
        cb.chk.show();
        cb.cover.isChecked = true;
    }
    cb.refresh = function(x) {
        if (x) {
            cb.check();
        }
        else {
            cb.uncheck();
        }
    }
    cb.hideIt = function() {
        cb.cover.hide();
        cb.box.hide();
        cb.chk.hide();
        for (var i = 0; i < cb.label.length; i++) {
            cb.label[i].hide();
        }
    }
    cb.showIt = function() {
        cb.cover.show();
        cb.box.show();
        if (cb.isChecked) {
            cb.chk.show();
        }
        else {
            cb.chk.hide();
        }
        for (var i = 0; i < cb.label.length; i++) {
            cb.label[i].show();
        }
    }
    var cbset = this.set(cb.box, cb.chk, cb.cover);
    cb.move = function(x, y) {
        cbset.transform("t" + x + "," + y);
    }
    return(cb);
}
Raphael.fn.radio = function(options) {
    if (missing(options.size)) {
        options.size = 6.5;
    }
    if (missing(options.vertspace)) {
        options.vertspace = 25;
    }
    if (missing(options.horspace)) { 
        options.horspace = rep(0, options.labels.length);
    }
    if (missing(options.lbspace)) {
        options.lbspace = 14;
    }
    if (missing(options.fontsize)) {
        options.fontsize = 14;
    }
    var rd = new Array();
    rd.whichChecked = options.whichChecked;
    rd.label = new Array(options.labels.length);
    rd.cover = new Array(options.labels.length);
    rd.circle = new Array(options.labels.length);
    rd.fill = this.set();
    var newvert = 0;
    for (var i = 0; i < options.labels.length; i++) {
        rd.label[i] = this.text(options.x + options.horspace[i] + options.lbspace, options.y + newvert - 1, options.labels[i]).attr({"text-anchor": "start", "font-size": options.fontsize+"px"});
        rd.circle[i] = this.circle(options.x + options.horspace[i], options.y + newvert, options.size).attr({fill: "#eeeeee", "stroke": "#a0a0a0", "stroke-width": 1.2});
        rd.cover[i] = this.circle(options.x + options.horspace[i], options.y + newvert, options.size + 2).attr({fill: "#eeeeee", stroke: "none", "fill-opacity": 0, "cursor": "pointer"});
        rd.cover[i].i = i;
        rd.cover[i].click(function() {
            rd.fill.show();
            var BBox = this.getBBox();
            rd.fill.transform("t" + (BBox.x - options.x + options.size + 2) + "," + (BBox.y - options.y + options.size + 2));
            rd.whichChecked = this.i;
        });
        if (Array.isArray(options.vertspace)) {
            newvert = options.vertspace[i + 1];
        }
        else {
            newvert = (i + 1)*options.vertspace;
        }
    }
    rd.fill.push(this.circle(options.x, options.y, options.size - 0.5).attr({fill: "#97bd6c", stroke: "none"}));
    rd.fill.push(this.circle(options.x, options.y, options.size - 4.5).attr({fill: "#000000", stroke: "none"}));
    if (rd.whichChecked < 0) {
        rd.fill.hide();
    }
    else {
        rd.fill.show();
        rd.fill.transform("t" + (rd.circle[options.whichChecked].getBBox().x - options.x + options.size) + "," + (rd.circle[options.whichChecked].getBBox().y - options.y + options.size));
    }
    rd.moveTo = function(pos) {
        rd.fill.show();
        rd.whichChecked = pos;
        var BBox = rd.cover[pos].getBBox();
        rd.fill.transform("t" + (BBox.x - options.x + options.size + 2) + "," + (BBox.y - options.y + options.size + 2));
    }
    rd.hideIt = function() {
        for (var i = 0; i < rd.cover.length; i++) {
            rd.cover[i].hide();
            rd.circle[i].hide();
            rd.label[i].hide();
        }
        rd.fill.hide();
    }
    rd.showIt = function() {
        for (var i = 0; i < rd.cover.length; i++) {
            rd.cover[i].show();
            rd.circle[i].show();
            rd.label[i].show();
        }
        if (rd.whichChecked < 0) {
            rd.fill.hide();
        }
        else {
            rd.fill.show();
        }
    }
    return(rd);
}
function missing(obj) {
    return(obj === void 0);
}
function isNumeric0(n) {
    if (n.length == 0) {
        return false;
    }
    else {
        return !/^(NaN|-?Infinity)$/.test(+n);
    }
}
function isNumeric(obj) {
    if (missing(obj)) {
        return false;
    } else if (obj === null) { 
        return false;
    } else if (obj.length == 0) {
        return false;
    }
    else {
        if (obj instanceof Array) {
            var result = new Array(obj.length);
            for (var i = 0; i < obj.length; i++) {
                result[i] = (obj[i].length == 0) ? false : !/^(NaN|-?Infinity)$/.test(+obj[i]);
            }
            return(result);
        }
        else {
            return !/^(NaN|-?Infinity)$/.test(+obj);
        }
    }
}
function copy(obj, exclude) {
    if (null == obj || "object" != typeof obj) return obj;
    var clone;
    if (obj instanceof Date) {
        copy = new Date();
        copy.setTime(obj.getTime());
        return copy;
    }
    if (obj instanceof Array) {
        clone = new Array();
        for (var i = 0, len = obj.length; i < len; i++) {
            if (exclude !== void 0) {
                if (exclude.indexOf(i) < 0) { 
                    clone.push(copy(obj[i]));
                }
            }
            else {
                clone.push(copy(obj[i]));
            }
        }
        return clone;
    }
    if (obj instanceof Object) {
        clone = new Object;
        for (var key in obj) {
            if (obj.hasOwnProperty(key)) {
                if (exclude !== void 0) {
                    if (exclude.indexOf(key) < 0) {
                        clone[key] = copy(obj[key]);
                    }
                }
                else {
                    clone[key] = copy(obj[key]);
                }
            }
        }
        return clone;
    }
}
function sortArray(obj, options) {
    if (missing(options)) {
        options = {};
    }
    if (options.emptylast == void 0) {
        options.emptylast = true;
    }
	var sorted = obj.slice();
	var sortlen = sorted.length;
	for (var i = 0, j; i < sortlen; i++) {
	    
	    tmp1 = sorted[i];
	    tmp2 = tmp1;
	    if (!isNumeric(tmp1)) {
	        tmp2 = tmp2.toLowerCase();
		}
		for (j = i - 1; j >= 0; j--) {
            tmp3 = sorted[j];
            if (!isNumeric(tmp3)) {
                tmp3 = tmp3.toLowerCase();
            }
            if (tmp3 > tmp2) {
                sorted[j + 1] = sorted[j];
            }
            else {
                break
            }
		}
		sorted[j + 1] = tmp1;
	}
	
	if (options.emptylast) {
	    obj = rep("", sortlen);
	    var position = 0;
	    for (var i = 0; i < sortlen; i++) {
	        if (sorted[i] != "") {
	            obj[position] = sorted[i];
	            position += 1;
	        }
	    }
	    return obj;
	}
	
	return sorted;
}
function objectsEqual(x, y) {
    'use strict';
    if (x === null || x === undefined || y === null || y === undefined) { return x === y; }
    if (x.constructor !== y.constructor) { return false; }
    if (x instanceof Function) { return x === y; }
    if (x instanceof RegExp) { return x === y; }
    if (x === y || x.valueOf() === y.valueOf()) { return true; }
    if (Array.isArray(x) && x.length !== y.length) { return false; }
    if (x instanceof Date) { return false; }
    if (!(x instanceof Object)) { return false; }
    if (!(y instanceof Object)) { return false; }
    var p = Object.keys(x);
    return Object.keys(y).every(function (i) { return p.indexOf(i) !== -1; }) &&
        p.every(function (i) { return objectsEqual(x[i], y[i]); });
}
function arraysEqual(a, b) {
    if (a === b) return true;
    if (a == null || b == null) return false;
    if (a.length != b.length) return false;
    for (var i = 0; i < a.length; ++i) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}
function getKeys(obj) {
    if (obj === null) return(Array());
    return(Object.keys(obj));
}
function getTrueKeys(obj) { 
    if (obj === null) {
        return(Array());
    }
    var trueKeys = new Array();
    for (var key in obj) {
        if (obj[key]) {
            trueKeys.push(key);
        }
    }
    if (trueKeys.length == 0) {
        return(Array());
    }
    else {
        return(trueKeys);
    }
}
function changeCol(obj, oldname, newname) {
    var temp = new Array();
    for (var key1 in obj) {
        temp[key1] = new Array();
        for (var key2 in obj[key1]) {
            temp[key1][key2] = new Array();
            var oldkeys = getKeys(obj[key1][key2]);
            for (var i = 0; i < oldkeys.length; i++) {
                temp[key1][key2][(oldkeys[i] == oldname)?(newname):(oldkeys[i])] = obj[key1][key2][oldkeys[i]];
            }
        }
    }
    return(temp);
}
function sat(obj, options) {
    if (options == void 0) {
        options = {};
    }
    if (options.size == void 0) {
        options.size = 14;
    }
    if (options.anchor == void 0) {
        options.anchor = "start";
    }
    var optkeys = getKeys(options);
    var options2 = new Array();
    for (var i = 0; i < optkeys.length; i++) {
        options2[optkeys[i]] = options[optkeys[i]];
    }
    if (Array.isArray(obj)) {
        var keys = getKeys(obj);
        for (var i = 0; i < keys.length; i++) {
            if (optkeys.indexOf("clip") >= 0) {
                if (typeof(options.clip) != "string") {
                    options2.clip = options.clip[keys[i]];
                }
            }
            sat(obj[keys[i]], options2);
        }
    }
    else {
        if (obj.type == "text") {
            var BBox, clipattr;
            if (options.clip == void 0) {
                obj.attr({"text-anchor": options.anchor, "font-size": (options.size + "px")});
                if (options["font-weight"] !== void 0) {
                    obj.attr({"font-weight": options["font-weight"]});
                }
                if (options["color"] !== void 0) {
                    obj.attr({"fill": options["color"]});
                }
            }
            else {
                if (typeof(options.clip) == "string") {
                    clipattr = options.clip;
                }
                else {
                    BBox = options.clip.getBBox();
                    clipattr = BBox.x + "," + BBox.y + "," + (BBox.width - 3) + "," + BBox.height;
                }
                obj.attr({"text-anchor": options.anchor, "font-size": (options.size + "px"), "clip-rect": clipattr});
            }
        }
        else { 
            obj.attr({stroke: ((optkeys.indexOf("stroke") < 0)?"#a0a0a0":options["stroke"]), 'stroke-width': ((optkeys.indexOf("sw") < 0)?1:options.sw), fill: "#ffffff", "fill-opacity": 0});
        }
    }
    return(obj);
}
function setAttrCover(obj) {
    if (Array.isArray(obj)) {
        var keys = getKeys(obj);
        for (var i = 0; i < keys.length; i++) {
            setAttrCover(obj[keys[i]]);
        }
    }
    else {
        obj.attr({stroke: '#d0d0d0', 'stroke-width': 1, fill: "#ffffff", "fill-opacity": 0})
    }
    return(obj);
}
function round(x, y) {
    y = Math.pow(10, y);
    return(Math.round(x*y)/y);
}
function all(obj, rule, value) {
    if (missing(value)) {
        value = "";
    }
    var check = true;
    if (Array.isArray(obj)) {
        var keys = getKeys(obj);
        for (var i = 0; i < obj.length; i++) {
            if (Array.isArray(value)) {
                var check2 = false;
                for (var j = 0; j < value.length; j++) {
                    check2 = check2 || eval("obj[keys[i]]" + rule + value[j]);
                }
                check = check && check2;
            }
            else {
                check = check && eval("obj[keys[i]]" + rule + value);
            }
        }
    }
    else {
        if (Array.isArray(value)) {
            var check = false;
            for (var j = 0; j < value.length; j++) {
                check = check || eval("obj" + rule + value[j]);
            }
        }
        else {
            check = eval("obj" + rule, value);
        }
    }
    return(check);
}
function any(obj, rule, value) {
    if (missing(value)) {
        value = "";
    }
    var check = false;
    if (Array.isArray(obj)) {
        var keys = getKeys(obj);
        for (var i = 0; i < obj.length; i++) {
            if (Array.isArray(value)) {
                for (var j = 0; j < value.length; j++) {
                    check = check || eval("obj[keys[i]]" + rule + value[j]);
                }
            }
            else {
                check = check || eval("obj[keys[i]]" + rule + value);
            }
        }
    }
    else {
        if (Array.isArray(value)) {
            for (var j = 0; j < value.length; j++) {
                check = check || eval("obj" + rule + value[j]);
            }
        }
        else {
            check = eval("obj" + rule + value);
        }
    }
    return(check);
}
function order(x, descending = false) {
    let values = new Array();
    for (let v = 0; v < x.length; v++) {
        values.push(x[v]);
    }
    for (var i = 0; i < values.length; i++) {
        values[i] = [values[i], i];
    }
    values.sort(function(left, right) {
        return((descending ? -1 : 1) * (left[0] < right[0] ? -1 : 1))
    });
    let result = [];
    for (var j = 0; j < values.length; j++) {
        result.push(values[j][1]);
    }
    return result;
}
function rep(rule, times) {
    var result = new Array(times);
    for (var i = 0; i < times; i++) {
        result[i] = rule;
    }
    return(result);
}
function unique(obj) {
    if (!Array.isArray(obj)) {
        return(null);
    }
    var uniques = new Array();
    var present;
    for (var i = 0; i < obj.length; i++) {
        if (uniques.indexOf(obj[i]) < 0) {
            uniques.push(obj[i]);
        }
    }
    return(uniques);
}
function min(obj) { 
    var minval = null;
    if (!Array.isArray(obj)) {
        return(minval);
    }
    for (var i = 0; i < obj.length; i++) {
        if (obj[i] !== null && isNumeric(obj[i])) {
            if (minval === null) {
                minval = obj[i];
            }
            else {
                if (minval > obj[i]) {
                    minval = obj[i];
                }
            }
        }
    }
    return(minval);
}
function max(obj) { 
    var maxval = null;
    if (!Array.isArray(obj)) {
        return(maxval);
    }
    for (var i = 0; i < obj.length; i++) {
        if (obj[i] !== null && isNumeric(obj[i])) {
            if (maxval === null) {
                maxval = obj[i];
            }
            else {
                if (maxval < obj[i]) {
                    maxval = obj[i];
                }
            }
        }
    }
    return(maxval);
}
function paste(obj, options) {
    if (Array.isArray(obj)) { 
        if (obj.length == 0) return("");
    }
    else {
        return("");
    }
    if (missing(options)) {
        options = {};
    }
    if (missing(options.sep)) {
        options.sep = " ";
    }
    if (missing(options.from)) {
        options.from = 0;
    }
    if (missing(options.to)) {
        options.to = obj.length - 1;
    }
    var result = obj[options.from];
    if (options.from < options.to) {
        for (var i = options.from + 1; i < options.to + 1; i++) {
            result += options.sep + obj[i];
        }
    }
    return(result);
}
function makeRules(oldv, newv) {
    var rule = new Array();
    if (oldv.length > 0) {
        rule = new Array(oldv.length);
        for (var i = 0; i < oldv.length; i++) {
            rule[i] = oldv[i] + "=" + newv[i];
        }
    }
    return(rule)
}
function eraseRecodeValues(paper) {
    paper.oldv.texts.VALUE.attr({"text": ""});
    paper.oldv.texts.range.FROM.attr({"text": ""});
    paper.oldv.texts.range.TO.attr({"text": ""});
    paper.oldv.texts.LOWESTTO.attr({"text": ""});
    paper.oldv.texts.TOHIGHEST.attr({"text": ""});
    paper.newv.texts.VALUE.attr({"text": ""});
}
function checkRecodeSelections(colclicks, paper) {
    var cols = getTrueKeys(colclicks.recode.rules);
    eraseRecodeValues(paper);
    if (cols.length == 1) {
        var lr = cols[0].split("=");
        var lhs = lr[0].split(":");
        var rhs = lr[1];
        var idx = 0; 
        if (lhs.length > 1) {
            if (lhs[0] == "lo") {
                idx = 1;
                paper.oldv.texts.LOWESTTO.attr({"text": lhs[1]});
                paper.rules.oldv = lhs[1];
            }
            else if (lhs[1] == "hi") {
                idx = 3;
                paper.oldv.texts.TOHIGHEST.attr({"text": lhs[0]});
                paper.rules.oldv = lhs[0];
            }
            else {
                idx = 2;
                paper.oldv.texts.range.FROM.attr({"text": lhs[0]});
                paper.oldv.texts.range.TO.attr({"text": lhs[1]});
                paper.rules.oldv = [lhs[0], lhs[1]];
            }
        }
        else {
            if (lhs[0] == "missing") {
                idx = 4;
            }
            else if (lhs[0] == "else") {
                idx = 5;
            }
        }
        if (idx == 0) { 
            paper.oldv.texts.VALUE.attr({"text": lhs[0]});
            paper.rules.oldv = lhs[0];
        }
        if (rhs != "missing" && rhs != "copy") {
            paper.newv.texts.VALUE.attr({"text": rhs});
            paper.newradio.moveTo(0);
            paper.rules.newv = rhs;
        }
        else {
            if (rhs == "missing") {
                paper.newradio.moveTo(1);
            }
            else {
                paper.newradio.moveTo(2);
            }
        }
        paper.oldradio.moveTo(idx);
        paper.add.attr({"text": "Change"});
    }
    else {
        paper.add.attr({"text": "Add"});
    }
}
function deleteRule(colclicks, rule) {
    var keys = getKeys(colclicks.recode.rules);
    var temp = new Array();
    for (i = 0; i < keys.length; i++) {
        if (keys[i] != rule) {
            temp[keys[i]] = colclicks.recode.rules[keys[i]];
        }
    }
    colclicks.recode.rules = temp;
}
function unselect(colclicks, dialog, identifier) {
    if (colclicks[dialog][identifier] != void 0) {
        var keys = getKeys(colclicks[dialog][identifier]);
        for (i = 0; i < keys.length; i++) {
            colclicks[dialog][identifier][keys[i]] = false;
        }
    }
}
function scaleplot(paper) {
    var xyplotdata = paper.xyplotdata;
    var scale = paper.scale;
    var sx = paper.sx;
    var sy = paper.sy;
    var dim = paper.dim;
    var offset = paper.offset;
    var rdim = paper.rdim;
    var xcoord, ycoord;
    var rj = paper.randomjitter;
    if (paper.total !== void 0) {
        paper.total.remove();
        paper.mdlines.remove();
        paper.afv.remove();
        paper.ticks.remove();
        paper.pointsset.remove();
    }
    paper.total = paper.set();
    paper.afv = paper.set(); 
    paper.ticks = paper.set();
    paper.mdlines = paper.set(); 
    paper.pointsset = paper.set();
    paper.points = new Array();
    paper.total.push(paper.rect(sx, sy, scale*dim, scale*dim));
    paper.total.push(paper.path("M" + sx + "," + (sy + scale*dim) + " L" + (sx + scale*dim) + "," + sy).attr({"stroke": "#a0a0a0"}));
    paper.mdlines.push(paper.path("M" + sx + "," + (sy + scale*dim/2) +  " L" + (sx + scale*dim) + "," + (sy + scale*dim/2)).attr({"stroke-dasharray": "--"}));
    paper.mdlines.push(paper.path("M" + (sx + scale*dim/2) + "," + sy +  " L" + (sx + scale*dim/2) + "," + (sy + scale*dim)).attr({"stroke-dasharray": "--"}));
    paper.mdlines.attr({"stroke": "#a0a0a0"});
    for (var i = 0; i < 11; i++) {
        var vertick = 
        paper.ticks.push(paper.path("M" + (sx + scale*(offset + i*rdim/10)) + "," + (sy + scale*dim) + " L" + (sx + scale*(offset + i*rdim/10)) + "," + (sy + scale*dim + 7)));
        paper.ticks.push(paper.path("M" + (sx - 7) + "," + (sy + scale*(offset + i*rdim/10)) + " L" + sx + "," + (sy + scale*(offset + i*rdim/10)) ) );
        paper.afv.push(sat(paper.text(sx + scale*(offset + i*rdim/10), sy + scale*dim + 15, i/10), {"size": 12, "anchor": "middle"}));
        paper.afv.push(sat(paper.text(sx - 10, sy + scale*(offset + i*rdim/10), (10 - i)/10), {"size": 12, "anchor": "end"}));
        if (i == 5) {
            paper.afv.push(sat(paper.text(sx + scale*(offset + i*rdim/10), sy + scale*dim + 34, ((paper.negx.isChecked?"~":"") + paper.x)), {"size": 14, "anchor": "middle", "font-weight": "bold"}));
            var temp = sat(paper.text(sx - 38, sy + scale*(offset + i*rdim/10), ((paper.negy.isChecked?"~":"") + paper.y)), {"size": 14, "anchor": "end", "font-weight": "bold"});
            var BBox = temp.getBBox();
            temp.transform("t" + (BBox.width/2) + ",0r-90");
            paper.afv.push(temp);
        }
    }
    var point, txt, txtfundal, totalength;
    var hoverIn = function() {
        if (this.label.attr("text") != "" & !paper.labels.isChecked) {
            this.label[1].attr({"fill-opacity": 0.4});
            this.label.show();
        }
    }
    var hoverOut = function() {
        if (this.label.attr("text") != "" & !paper.labels.isChecked) {
            this.label[1].attr({"fill-opacity": 0});
            this.label.hide();
        }
    }
    if (xyplotdata.length > 0) { 
        var krj = Object.keys(rj);
        var px = xyplotdata[1];
        var py = xyplotdata[2];
        for (var i = 0; i < xyplotdata[0].length; i++) {
            xcoord = sx + scale*(offset + rdim*((paper.negx.isChecked)?(1 - px[i]):(px[i]))) + ((krj.length > 0)?(rj.x[i]):0);
            ycoord = sy + scale*(offset + rdim*((paper.negy.isChecked)?(py[i]):(1 - py[i]))) + ((krj.length > 0)?(rj.y[i]):0);
            paper.points[i] = paper.circle(xcoord, ycoord, 3);
            paper.points[i].hover(hoverIn, hoverOut, point, point);
            paper.pointsset.push(paper.points[i]);
        }
        if (paper.pof.isChecked) {
            if (paper.sufnec.whichChecked == 0) {
                paper.incl.attr({"text": "Inclusion: " + xyplotdata[3][paper.index][0]});
                paper.cov.attr({"text": "Coverage: " + xyplotdata[3][paper.index][1]});
                paper.PRI.attr({"text": "PRI: " + xyplotdata[3][paper.index][2]});
                paper.ron.hide();
            }
            else {
                paper.incl.attr({"text": "Inclusion: " + xyplotdata[4][paper.index][0]});              
                paper.cov.attr({"text": "Coverage: " + xyplotdata[4][paper.index][1]});
                paper.ron.attr({"text": "Relevance: " + xyplotdata[4][paper.index][2]});
                paper.PRI.hide()
            }
        }
        paper.pointsset.attr({fill: "#707070", "fill-opacity": paper.fill.isChecked?1:0});
    }
}
function createLabels(paper) {
    var xyplotdata = paper.xyplotdata;
    paper.labelsset.remove();
    paper.labelsset = paper.set();
    paper.labelsArray = new Array();
    var scale = paper.scale;
    var sx = paper.sx;
    var sy = paper.sy;
    var dim = paper.dim;
    var offset = paper.offset;
    var rdim = paper.rdim;
    var txt, txtfundal, outer, BBox, coords, x, y, r, twidth;
    if (xyplotdata.length > 0) { 
        for (var i = 0; i < xyplotdata[0].length; i++) {
            twidth = getTextWidth(xyplotdata[0][i])
            x = paper.points[i].attr("cx");
            y = paper.points[i].attr("cy");
            r = paper.points[i].attr("r");
            outer = paper.circle(x, y, r + twidth + 2*5).attr({"fill": "none", "stroke": "none"});
            coords = paper.points[i].getPointAtLength(paper.points[i].getTotalLength()*(90 - paper.labelRotation)/360);
            txtfundal = paper.rect(coords.x, coords.y - 8, twidth + 10, 16);
            txtfundal.attr({fill: "#c9c9c9", "fill-opacity": 0, stroke: "none"});
            txt = sat(paper.text(coords.x + 2.5, coords.y, xyplotdata[0][i]));
            txt.attr({"font-weight": "bold", "fill-opacity": 0.7});
            tempset = paper.set(paper.points[i], txt, txtfundal);
            BBox = tempset.getBBox();
            if (BBox.x + BBox.width > sx + scale*rdim) {
                txt.remove();
                txtfundal.remove();
                coords = outer.getPointAtLength(outer.getTotalLength()*(270 - paper.labelRotation)/360);
                txtfundal = paper.rect(coords.x, coords.y - 8, twidth + 10, 16);
                txtfundal.attr({fill: "#c9c9c9", "fill-opacity": 0, stroke: "none"});
                txt = sat(paper.text(coords.x + 2.5, coords.y, xyplotdata[0][i]));
                txt.attr({"font-weight": "bold", "fill-opacity": 0.7});
            }
            tempset = paper.set(paper.points[i], txt, txtfundal);
            paper.points[i].label = paper.set();
            paper.points[i].label.push(txt, txtfundal);
            if (paper.labelRotation > 0) {
                paper.points[i].label.transform("r-" + paper.labelRotation + "," + coords.x + "," + coords.y);
            }
            paper.labelsset.push(paper.points[i].label);
        }
    }
    if (paper.labels.isChecked) {
        paper.labelsset.show();
    }
    else {
        paper.labelsset.hide();
    }
}
function setPath(paper, wd, width) {
    paper.clear();
    var wds = wd.split("/");
    var x = 1;
    var dirwidth = 0;
    var direct;
    var allobjs = paper.set();
    for (i = 0; i < wds.length; i++) {
        if (wds[0] == "") {
            wds[0] = "root";
        }
        if (wds[i] != "") {
            dirwidth = getTextWidth(wds[i]);
            allobjs.push(sat(paper.text(x, 7, wds[i])));
            direct = paper.rect(x - 3, 0, dirwidth + 6, 18)
                .attr({stroke: "none", fill: "#fff", "fill-opacity": 0});
            direct.dir = wds[i];
            direct.dblclick(function() {
                paper.goToDir(this.dir)
            });
            allobjs.push(direct);
            x += dirwidth + 2;
            if (i < wds.length - 1*((wds[wds.length - 1] == "")?2:1)) {
                allobjs.push(sat(paper.text(x + 3 , 7, ">")));
                x += 15;
            }
        }
    }
    var lastX, groupX, groupWidth;
    function dragStartLarge(group) {
        return function() {
            lastX = 0;
            var BBox = group.getBBox();
            groupX = BBox.x;
            groupWidth = BBox.width;
        }
    };
    function dragMoveLarge(group) {
        return function(dx, dy) {
            if (groupWidth > (width + 50)) {
                var newX = dx - lastX;
                if (groupX + dx > 0) {
                    newX = 0 - groupX - lastX;
                }
                if (groupX + dx < (width - groupWidth)) {
                    newX = (width - groupWidth) - (groupX + lastX);
                    lastX = (width - groupWidth) - groupX;
                }
                else {
                    lastX += newX;
                }
                group.translate(newX, 0);
            }
        }
    };
    function dragStopLarge(group) {
        return function() {
        }
    };
    allobjs.drag(dragMoveLarge(allobjs), dragStartLarge(allobjs), dragStopLarge(allobjs));
    paper.setSize(x + 20, 20);
}
function randomBetween(min, max) {
    return (min + Math.random()*(max - min + 1));
}
function reorder(obj, from, to) {
    var keys = getKeys(obj);
    var values = new Array(keys.length);
    for (var i = 0; i < keys.length; i++) {
        values[i] = obj[keys[i]];
    }
    keys.splice(to, 0, keys.splice(from, 1)[0]);
    values.splice(to, 0, values.splice(from, 1)[0]);
    var result = {};
    for (var i = 0; i < keys.length; i++) {
        result[keys[i]] = values[i];
    }
    return(result);
}
function duplicates(arr) {
    var len = arr.length,
        out = [],
        counts = {};
    for (var i = 0; i < len; i++) {
        var item = arr[i];
        counts[item] = counts[item] >= 1 ? counts[item] + 1 : 1;
    }
    for (var item in counts) {
        if (counts[item] > 1) {
            out[out.length] = item*1;
        }
    }
    return out;
}
function decToBin(x) {
    if (isNumeric(x)) {
        return((_$=($,_="")=>$?_$($>>1,($&1)+_):_)(x));
    } else {
        return("");
    }
}
function getScrollBarWidth() {
    var inner = document.createElement('p');
    inner.style.width = "100%";
    inner.style.height = "200px";
    var outer = document.createElement('div');
    outer.style.position = "absolute";
    outer.style.top = "0px";
    outer.style.left = "0px";
    outer.style.visibility = "hidden";
    outer.style.width = "200px";
    outer.style.height = "150px";
    outer.style.overflow = "hidden";
    outer.appendChild (inner);
    document.body.appendChild (outer);
    var w1 = inner.offsetWidth;
    outer.style.overflow = 'scroll';
    var w2 = inner.offsetWidth;
    if (w1 == w2) w2 = outer.clientWidth;
    document.body.removeChild (outer);
    return (w1 - w2);
};
function scaleShape(path, scale) {
    var parsed = Raphael.parsePathString(path);
    for (var j = 0; j < parsed.length; j++) {
        for (k = 1; k < parsed[j].length; k++) {
            parsed[j][k] = parsed[j][k]*scale;
        }
    }
    return(parsed.toString());
}
function getShape(x, venn, scale) {
    bigpath = "";
    for (var b = 0; b < x.length; b++) {
        var path = "M";
        var stb, endb, end;
        var checkb = rep(false, x[b].length);
        var counter = 0;
        while(counter < 1000) { 
            for (var i = 0; i < checkb.length; i++) {
                if (!checkb[i]) {
                    var y = venn[x[b][i]];
                    if (i == 0) {
                        for (var j = 0; j < y.length/2; j++) {
                            path += ((j == 1)?" C":"") + " " + round(y[2*j]*scale, 3) + "," + round(y[2*j + 1]*scale, 3);
                        }
                        checkb[i] = true;
                        end = y[y.length - 2] + " " + y[y.length - 1];
                    }
                    else {
                        stb = y[0] + " " + y[1];
                        endb = y[y.length - 2] + " " + y[y.length - 1]
                        if (end == stb) {
                            for (var j = 1; j < y.length/2; j++) {
                                path += " " + round(y[2*j]*scale, 3) + "," + round(y[2*j + 1]*scale, 3);
                            }
                            checkb[i] = true;
                            end = endb;
                        }
                        else if (end == endb) {
                            for (var j = y.length/2 - 2; j >= 0; j--) {
                                path += " " + round(y[2*j]*scale, 3) + "," + round(y[2*j + 1]*scale, 3);
                            }
                            checkb[i] = true;
                            end = stb;
                        }
                    }
                }
            }
            if (all(checkb, "== true")) {
                counter = 1001;
            }
            counter += 1
        }
        bigpath += " " + path + " z";
    }
    return(bigpath);
}
function customShape(rule, venn, scale, id) {
    var rowns = new Array();
    rule = rule.split("");
    var idis, i, j, k;
    var check = rep(true, rule.length);
    var keys = getKeys(id);
    for (i = 0; i < keys.length; i++) { 
        keys[i] = keys[i]*1;
        idis = id[keys[i]].split("");
        for (j = 0; j < rule.length; j++) {
            if (rule[j] != "-") {
                check[j] = any(idis, "==", j + 1);
                if (rule[j] == "0") {
                    check[j] = !check[j];
                }
            }
        }
        if (all(check, "== true")) {
            rowns[rowns.length] = i;
        }
    }
    var ids = new Array();
    for (i = 0; i < rowns.length; i++) {
        ids[i] = id[rowns[i]];
    }
    var inverted = any(rowns, "==", 0);
    if (rowns.length == 1 && rowns[0] == 0) {
        rowns = $(keys).not(rowns).get();
    }
    checkZone = function(from, rowns, checkz, venn) {
        var fromz = venn[1][from];
        var toz = new Array();
        for (var i = 0; i < rowns.length; i++) {
            if (!checkz[i]) {
                if (any(fromz, "==", venn[1][rowns[i]])) {
                    checkz[i] = true;
                    toz[toz.length] = rowns[i];
                }
            }
        }
        if (toz.length > 0) {
            for (var j = 0; j < toz.length; j++) {
                var checkz2 = checkZone(toz[j], rowns, checkz, venn);
                for (var i = 0; i < checkz.length; i++) {
                    checkz[i] = checkz[i] || checkz2[i];
                }
            }
        }
        return(checkz)
    }
    var result = new Array();
    if (rowns.length > 1) {
        var checkz = rep(false, rowns.length);
        checkz[0] = true;
        while(any(checkz, "== false")) {
            checkz = checkZone(rowns[0], rowns, checkz, venn);
            var temp1 = new Array();
            var temp2 = new Array();
            var checkz2 = new Array();
            for (i = 0; i < rowns.length; i++) {
                if (checkz[i]) {
                    temp1[temp1.length] = rowns[i];
                }
                else {
                    temp2[temp2.length] = rowns[i];
                    checkz2[checkz2.length] = false;
                }
            }
            result[result.length] = temp1;
            if (checkz2.length > 0) {
                rowns = copy(temp2);
                checkz = copy(checkz2);
                checkz[0] = true;
            }
        }
    }
    else {
        result[0] = [rowns];
    }
    for (var i = 0; i < result.length; i++) {
        var temp = venn[1][result[i][0]];
        if (result[i].length > 1) {
            for (var j = 1; j < result[i].length; j++) {
                temp = temp.concat(venn[1][result[i][j]]);
            }
        }
        result[i] = $(temp).not(duplicates(temp)).get();
    }
    return([getShape(result, venn[0], scale), inverted]);
}
function parseText(text, conditions) {
    text = text.replace("(", "");
    text = text.replace(")", "");
    text = text.replace(/\s/g, "");
    splitchar = "*";
    var parsedPlus = text.split("+");
    var largecheck = rep(false, parsedPlus.length);
    for (var i = 0; i < parsedPlus.length; i++) {
        var parsedStar = parsedPlus[i].split(splitchar);
        var upper = new Array(parsedStar.length);
        var lower = new Array(parsedStar.length);
        for (var j = 0; j < parsedStar.length; j++) {
            if (parsedStar[j][0] == "~") {
                parsedStar[j] = parsedStar[j].substring(1, parsedStar[j].length);
            }
            upper[j] = parsedStar[j].toUpperCase();
            lower[j] = parsedStar[j].toLowerCase();
            if (parsedStar[j] != upper[j] && parsedStar[j] != lower[j]) {
                return("error");
            }
        }
        var check = rep(false, parsedStar.length);
        for (var j = 0; j < parsedStar.length; j++) {
            check[j] = conditions.indexOf(upper[j]) < 0;
        }
        largecheck[i] = any(check, " == true")
    }
    if (any(largecheck, " == true")) {
        splitchar = "";
        var largecheck = rep(false, parsedPlus.length);
        for (var i = 0; i < parsedPlus.length; i++) {
            var parsedStar = parsedPlus[i].split(splitchar);
            var upper = new Array(parsedStar.length);
            for (var j = 0; j < parsedStar.length; j++) {
                if (parsedStar[j][0] == "~") {
                    parsedStar[j] = parsedStar[j].substring(1, parsedStar[j].length);
                }
                upper[j] = parsedStar[j].toUpperCase();
            }
            var check = rep(false, parsedStar.length);
            for (var j = 0; j < parsedStar.length; j++) {
                check[j] = conditions.indexOf(upper[j]) < 0;
            }
            largecheck[i] = any(check, " == true")
        }
        if (any(largecheck, " == true")) {
            return("error");
        }
    }
    finalResult = {};
    for (var i = 0; i < parsedPlus.length; i++) {
        var parsedStar = parsedPlus[i].split(splitchar);
        var upper = new Array(parsedStar.length);
        for (var j = 0; j < parsedStar.length; j++) {
            if (parsedStar[j][0] == "~") {
                parsedStar[j] = parsedStar[j].substring(1, parsedStar[j].length).toLowerCase();
            }
            upper[j] = parsedStar[j].toUpperCase();
        }
        var rule = "";
        for (var j = 0; j < conditions.length; j++) {
            var index = upper.indexOf(conditions[j]);
            rule += (index >= 0)?((parsedStar[index] == upper[index])?1:0):"-";
        }
        finalResult[parsedPlus[i]] = rule;
    }
    return(finalResult);
}
function caretPosition(jqueryItem) {
    var input = jqueryItem.get(0);
    if (!input) return; 
    if ('selectionStart' in input) {
        return input.selectionStart;
    } else if (document.selection) {
        input.focus();
        var sel = document.selection.createRange();
        var selLen = document.selection.createRange().text.length;
        sel.moveStart('character', -input.value.length);
        $("#sel").html(sel);
        $("#selLen").html(selLen);
        return sel.text.length - selLen;
    }
}
function getStyle(fresh) {
    if (missing(fresh)) {
        fresh = false;
    }
    var oStyle = {
        background: 'none',
        width: (fresh?651:(($("#txtcommand").length)?($("#txtcommand").width()):651)) + 'px',
        height: (fresh?19:(($("#txtcommand").length)?($("#txtcommand").height()):19)) + 'px',
        zIndex: '9000',
        padding: '0 0 0 0',
        border: 'none', 
        resize: 'none',
        outline: 'none',
        'font-size': '14px',
        'font-family': "Monaco,Menlo,Consolas,'Courier New',monospace",
        color: 'blue',
        'overflow-y': 'hidden'
    }
    var sStyle = '';
    for (var z in oStyle){
        sStyle += z + ':' + oStyle[z] + ';';
    }
    return(sStyle);
}
function escapeText( s ) {
	if ( !s ) {
		return "";
	}
	s = s + "";
	
	return s.replace( /['"<>&]/g, function( s ) {
		switch ( s ) {
		case "'":
			return "&#039;";
		case "\"":
			return "&quot;";
		case "<":
			return "&lt;";
		case ">":
			return "&gt;";
		case "&":
			return "&amp;";
		}
	});
}
