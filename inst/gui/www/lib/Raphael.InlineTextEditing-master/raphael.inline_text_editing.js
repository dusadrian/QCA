/*
 * Inline text editing tool for Raphaël 2.0 & compatible with Raphaël Free transform.
 * Source: https://github.com/marmelab/Raphael.InlineTextEditing
 * Licensed under the MIT license
 * Modified by Adrian Dusa, 2015
 */

(function (root, factory) {
	if (typeof define === "function" && define.amd) {
		// AMD. Register as an anonymous module.
		define(["raphael"], function(Raphael) {
			// Use global variables if the locals are undefined.
			return factory(Raphael || root.Raphael);
		});
	} else {
		// RequireJS isn't being used. Assume Raphael is loaded in <script> tag
		factory(Raphael);
	}
}(this, function(Raphael) {
    
    var typeofcall;
    

	Raphael.fn.inlineTextEditing = function(subject, options, callback) {
	    
	    // Store instance of the Raphael paper
		var paper = this;

		subject.inlineTextEditing = {
			paper : paper,
			input: null,

			/**
			 * Start text editing by hiding the current element and adding a text field at the same position
			 * @return jQuery input element
			 */
			startEditing: function(x, y, width, height, type, background) {
			    /* AD
			    the argument "type" was added because only the data editor is supposed to
			    change the width of the textbox
			    */
			    typeofcall = type;
			    if (background === undefined) {
			        background = 'white';
			    }
			    
				// Store Raphael container above the svg
				// modification AD (2015.11.10):
				// instead of .parentNode is now .parentNode.parentNode
				// to point to the outer dialog, not to the subdiv of that dialog
				var container       = this.paper.canvas.parentNode.parentNode;
				var translateX	   = 0;
				var translateY	   = 0;
				var transformOrder  = {};

				// Retrieve element transformation
				var rotation        = subject._.deg;
				var scaleX          = subject._.sx;
				var scaleY          = subject._.sy;
				var matrix          = subject.node.getAttribute('transform');

				// Check if the element has translations & retrieve transformations order
				for(var i = 0, length = subject._.transform.length; i < length; i++){
					var matrixComponents = subject._.transform[i];
					var transform = matrixComponents[0].toLowerCase();
					transformOrder[transform] = transform;

					if(transform == 't'){
						translateX += matrixComponents[1];
						translateY += matrixComponents[2];
					}
				}


				// Check if there is implicit matrix
				for(var i = 0, length = subject._.transform.length; i < length; i++){
					if(subject._.transform[i][0].toLowerCase() == 'm'){
						var matrixComponents = subject._.transform[i].slice(1);

						// Perform transformation from matrix elements
						rotation  += -1 * Math.asin(matrixComponents[2]) * 180 / Math.PI;
						scaleX    *= matrixComponents[0] / Math.cos(rotation*Math.PI/180);
						scaleY    *= matrixComponents[3] / Math.cos(rotation*Math.PI/180);

						transformOrder = {r:'r', s:'s'};
					}
				}

				// Remove transformation on the current element to retrieve original dimension
				subject.node.removeAttribute('transform');

				//var originalBbox  = subject._getBBox();
				//var width         = 70;
				//var height        = 20;
				//var x             = originalBbox.x;
				//var y             = originalBbox.y;
				var sTransform    = '';
				var sOrigin       = 'center center';
				var oTransform    = {
					//	t : 'translate('+(translateX)+'px, '+(translateY)+'px)',
					r : 'rotate('+rotation+'deg)',
					s : 'scale('+scaleX+', '+scaleY+')'
				};

				// Build transform CSS property in the same order than the element
				for (var transform in transformOrder){
					if(oTransform[transform] != undefined){
						sTransform += oTransform[transform] + ' ';
					}
				}

				// Re-apply stored transformation to the element and hide it
				//subject.node.setAttribute("transform", matrix);
				subject.hide();
				
				// Prepare input styles
				var oStyles = {
					position: 'absolute',
					background: background,
					left: x+'px',
					top: y+'px',
					width: width+'px',
					height: height+'px',
					color: subject.attrs.fill,
					zIndex: '9000',
					padding: ((background == 'white')?'1 0 0 4':'2 0 0 5'),
					border: ((background == 'white')?'none':'#d7d7d7'),
					resize: 'none',
					outline: 'none',
					'overflow-y': 'hidden'
					/*
					'-moz-transform-origin': sOrigin,
					'-ms-transform-origin': sOrigin,
					'-o-transform-origin': sOrigin,
					'-webkit-transform-origin': sOrigin,
					'transform-origin': sOrigin,

					'-moz-transform' : sTransform,
					'-ms-transform' : sTransform,
					'-o-transform' : sTransform,
					'-webkit-transform' : sTransform,
					'transform' : sTransform
					*/
				};

				// Retrieve font styles
				var aFontAttributes = ['font', 'font-family', 'font-size', 'font-style', 'font-weight', 'font-variant'/*, 'line-height'*/];

				for (var i = 0; i < aFontAttributes.length; i++){
					var attribute = aFontAttributes[i];

					if (subject.attrs[attribute] != undefined){
						oStyles[attribute] = subject.attrs[attribute];
					}

					if (subject.node.style[attribute] != undefined){
						oStyles[attribute] = subject.node.style[attribute];
					}
				}

				var sStyles = '';
				for (var z in oStyles){
					sStyles += z + ':' + oStyles[z] + ';';
				}
				
				
				
				// Create an input element with theses styles
				this.input = document.createElement("textarea");
				//contenteditable la div face o chestie
				
				var txt = "";
				if (subject.attrs.text) {
				    txt = subject.attrs.text.replace(/\'/g,"\\\'");
				}
                 
                 // if (type == "from_filldirexp") {
                 //    if (txt == "-") {
                 //        txt = "";
                 //    }
                 // }
                 
                 this.input.value = txt;
                 
                 this.input.setAttribute("style", sStyles);
                 
				this.input.addEventListener('keyup', this._handleKeyDown.bind(this));

				this.input.id = "txtarea";
				
				// Add the input in the container and apply focus on it
				container.appendChild(this.input);
				
				
				$("#txtarea").click(function(event) {
                    event.stopPropagation();
				});
				
				this.input.focus();
				
				
				return this.input;
			},

			/**
			 * Apply text modification and remove associated input
			 */
			stopEditing: function(key) {

				// Set the new value
				if (key == "enter") {
				    subject.attr("text", this.input.value);
				}
				
				// Show the text element
				subject.show();

				// Remove text input
				this.input.parentNode.removeChild(this.input);
			},
			
			_handleKeyDown: function(e){
				var tmp               = document.createElement("span");
				var text              = this.input.value;
				tmp.setAttribute('style', this.input.style.cssText);
				tmp.style.height      = null;
				tmp.style.width       = null;
				tmp.style.visibility  = 'hidden';
				//tmp.innerHTML         = text.split('\n').join('<br />');

				this.input.parentNode.appendChild(tmp);
				
				var textwidth = getTextWidth(text);
				/*
				// only for the data editor does the textbox width changes
				// 
				if (typeofcall == "from_data_editor") {
				
                     if (textwidth >= 40) {
                         textwidth += 30
                     }
                     else {
                         textwidth = 70;
                     }
                     
                     this.input.style.width = textwidth + "px";
                     //this.input.style.height = (tmp.offsetHeight + 10) + "px";this.input.addEventListener("keypress", function() {
                    
                    
                 }
                 */
                 if (this.input.scrollTop != 0) {
                    this.input.style.height = this.input.scrollHeight + "px";
                 }

				tmp.parentNode.removeChild(tmp);
			}
		};

		return subject.inlineTextEditing;
	}

}));


(function( $ ) {
	$.fn.stopScroll = function( options ) {
		options = $.extend({
			delay: 250,
			callback: function() {}
		}, options);
		
		return this.each(function() {
			var $element = $( this ),
				element = this;
			$element.scroll(function() {
				clearTimeout( $.data( element, "scrollCheck" ) );
				$.data( element, "scrollCheck", setTimeout(function() {
					options.callback();
				}, options.delay ) );
			});
		});
	};

})( jQuery );



/*
(function( $ ) {
	$.fn.Rcommand = function(subject, options, callback) {
	    console.log(subject.selector);
	    subject.Rcommand = {
            
            input: null,
            
            startEditing: function(position) {
                var container = this.parentNode;
                var oStyles = {
                    position: 'absolute',
                    background: 'none',
                    left: (position.left + 15) + 'px',
                    top: position.top + 'px',
                    width: ($("#" + divid).width() - 15) + 'px',
                    height: $("#" + divid).height() + 'px',
                    zIndex: '9000',
                    padding: '1 0 0 4',
                    border: '1', // 'none'
                    resize: 'none',
                    outline: 'none',
                    'font-size': '14px',
                    'font-family': "Monaco,Menlo,Consolas,'Courier New',monospace"
                }
                
                
                var sStyles = '';
                for (var z in oStyles){
                    sStyles += z + ':' + oStyles[z] + ';';
                }
                
                this.input = document.createElement("textarea");
                this.input.value = "";
                
                this.input.setAttribute("style", sStyles);
                this.input.addEventListener('keyup', this._handleKeyDown.bind(this));
                
                this.input.id = "txtarea";
                
                container.appendChild(this.input);
                
                
				$("#txtarea").click(function(event) {
                    event.stopPropagation();
				});
				
				this.input.focus();
				
				return this.input;
            },
            stopEditing: function(key) {
                
                // Set the new value
                if (key == "enter") {
                    subject.attr("text", this.input.value);
                }
                
                // Remove text input
                this.input.parentNode.removeChild(this.input);
            },
            _handleKeyDown: function(e) {
                var tmp               = document.createElement("span");
                var text              = this.input.value;
                tmp.setAttribute('style', this.input.style.cssText);
                tmp.style.height      = null;
                tmp.style.width       = null;
                tmp.style.visibility  = 'hidden';
                //tmp.innerHTML         = text.split('\n').join('<br />');
    
                this.input.parentNode.appendChild(tmp);
                
                
                var textwidth = getTextWidth(text);
                if (this.input.scrollTop != 0) {
                    this.input.style.height = this.input.scrollHeight + "px";
                }
    
                tmp.parentNode.removeChild(tmp);
            }
	    }
	    
	    return subject.editText;
	};

})( jQuery );
*/


function getTextWidth(string) {
    var paper = Raphael(0, 0, 0, 0);
        paper.canvas.style.visibility = "hidden";
    var el = paper.text(0, 0, string).attr({"text-anchor": "start", "font-size": "14px"});
    var bBox = el.getBBox();
    paper.remove();
    return bBox.width;
}




