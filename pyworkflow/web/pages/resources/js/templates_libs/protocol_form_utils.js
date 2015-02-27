 /*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
 * 			   Adrian Quintana (aquintana@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'jmdelarosa@cnb.csic.es'
 *
 ******************************************************************************/
/******************************************************************************
 * DESCRIPTION:
 * 
 * Methods to manage the protocol form.
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * 
 * jQuery(document).ready(function())
 * 	->	Overray the post simple method in the protocol form template. 
 * 		Depend a variable of the protocol form, can be:
 * 		* Execute protocol: This mode execute the protocol with the parameters of
 * 			the form.
 * 		* Save protocol: This method save the protocol with the parameters filled.
 * 		* Wizard: This mode launch a wizard for a specific parameter.   
 * 		* Viewer: This mode launch a viewer to analyze the results for the parameters
 * 			chosen in the form.
 * 
 * function fixInput(serialize_form)
 * 	->	Function to replace the label for a input object by his objId in the 
 * 		serialized form with his parameters.
 * 
 * function evalElements()
 * 	->	Function to evaluate the elements in a form, depending the type of the
 * 		param, it is evaluated in a diferent way.
 * 
 * function onChangeParam(value, paramId)
 * 	->	Update the parameter for an element.
 * 
 * function onChangeEnumParamCombo(elemId, paramId)
 * 	->	Update the parameter for an element type ENUM_PARAM_COMBO.
 * 
 * function onChangeEnumParamList(index, paramId)
 * 	->	Update the parameter for an element type ENUM_PARAM_LIST.
 * 
 * function setParamValue(paramId, value)
 * 	->	Put the new value in an attribute of the parent node.
 * 
 * function evalDependencies(row, newLevel)
 * 	->	Function to evaluate the parameters dependencies for a expertise level.
 * 
 * function evalCondition(row)
 * 	->	Function to evaluate a condition given a row of the form.
 * 
 * function normalizeConditions(cond)
 * 	->	For some conditions not normalized, it is replaced to be evaluated right.
 * 
 * function browseObjects(param, projName, type_param, param)
 * 	->	Browse object in the database depending of a type_param. 
 * 		Params: 
 * 			objClass: the class to get instances from (also subclasses)
 * 			protClassName: class refered to a protocol
 * 
 * function formProtSimple(param, projName)
 * 	->	Launch a custom protocol form with less options, thought for workflows
 * 		where some options not need to be chosen.
 * 
 * function returnProtocol()
 * 	->	Function to return a simple protocol form with the values serializes to be
 * 		saved, used for the custom workflow where the parameters need to be changed.
 * 
 * function setParamProt(paramProt, params)
 * 	->	Function to update a protocol serialize param inside the workflow.
 * 		The serializes values are storaged inside the input respectively.
 * 
 * function getTableFormatted(node, list, id, previsualize)
 * 	->	This function return an html table with the elements available after use, 
 * 		the function browse_object.
 * 
 * function selectDialog(objClass, msg, funcName)
 * 	->	This function create a popup with messi.js to choose the elements showed
 * 		with the method getTableFormatted.
 * 
 * function processSelectionTable(elm)
 * 	->	After select an element (PointerParam) showed with the function selectDialog with the 
 * 		table generated by the function getTableFormatted, this selection is 
 * 		processed with this method.
 * 
 * function processMultipleSelectionTable(elm)
 * 	->	After select an element (MultiPointerParam) showed with the function selectDialog with the 
 * 		table generated by the function getTableFormatted, this selection is 
 * 		processed with this method.
 * 
 * function removeObjects(paramName)
 *  ->  Remove element/s from multiple select element (it is usually associated to a MultiPointerParam)
 * 
 * function recalculateSelectDim(selectElementName)
 *  ->  Calculate dim for select multiple element depending on the number of options.
 * 		This is usually called after adding or removing an option
 * 
 * function selTableMessi(elm)
 * 	->	Function to select the element dinamically with flashlighted colors in 
 * 		the table generated with the function getTableFormatted.
 * 
 ******************************************************************************/

/** METHODS ******************************************************************/

$(document).ready(function() {
	/*	
	* Overray the post simple method in the protocol form template. 
	* Depend a variable of the protocol form 
	*/
	$("#protocolForm").submit(function() {
		var mode = $("#protocolForm").attr('data-mode');

// ---------- MODE EXECUTE PROTOCOL  -----------------------------------------------
		if (mode == 'execute') {
			/* Execute the protocol */
			// console.log($("#protocolForm").serialize())
			var action = $("#protocolForm").attr("action");
			var URL = getSubDomainURL() + action
			var serialize_form = fixInput($("#protocolForm").serialize());
			
			$.post(URL, serialize_form, function(json) {
				if (json.errors.length > 0) {
					// Show errors in the validation
					errorPopup('Errors found', json.errors);
				} else {
					infoPopup('Success', "The protocol was launched successfuly",1);
				}
			}, "json");

// ---------- MODE SAVE PROTOCOL ------------------------------------------------
		} else if (mode == 'save') {
			/* Save the protocol */
			var serialize_form = fixInput($("#protocolForm").serialize());
			var URL = getSubDomainURL() + "/save_protocol/"
			$.post(URL, serialize_form, function(json) {
				if (json.errors != undefined) {
					// Show errors in the process to save
					errorPopup('Errors found',json.errors);
				} else {
					// No errors in the process to save
					protId = json.success;
					infoPopup('Success', "The protocol was saved successfuly", 1, 'window.opener.popup(\'/form/?protocolId='+protId+'\')');
				}
			},"json");

// ---------- MODE WIZARD -----------------------------------------------------
		} else if (mode == 'wiz') {
			
			new Messi("<i class='fa fa-magic'/>  Loading Wizard...",{
				modal : true
				});
			
			/* Execute the wizard */
			var serialize_form = fixInput($("#protocolForm").serialize());
			var URL = getSubDomainURL() + "/wizard/"
			$.post(URL, serialize_form, function(html) {
				
				$('.messi').remove();
				$('.messi-modal').remove();
				
				if(html=="errorInput"){		// Important. Stop the normal POST
					errorPopup("Error", "Input was not selected, please choose one.");
				} else if (html=="errorEmpty"){
					errorPopup("Error", "Input selected are None");
				} else if (html=="errorIterate"){
					errorPopup("Error", "Error iterating over the input set");
				} else {
					customPopupHTML(html,800,540);
				}
				return false;
			});

// ---------- MODE VIEWER -----------------------------------------------------
		} else if (mode == 'viewerElement') {
			
			new Messi("<i class='fa fa-eye'/> Loading Viewer...",{
				modal : true
				});
			
			/* Launch the viewers with the options chosen */
			var serialize_form = fixInput($("#protocolForm").serialize());
			var URL = getSubDomainURL() + "/viewer_element/"
			$.post(URL, serialize_form, function(json) {
				$('.messi').remove();
				$('.messi-modal').remove();				
				popUpJSON(json);
			},"json");			
		} 
		// Important. Stop the normal POST
		return false;
	});
});

function fixInput(serialize_form){
	var attrs = serialize_form.split("&");

	$.each(attrs, function(param, paramName) {
		// console.log(param)
		// console.log(paramName)
		var aux = paramName.split("=");
		if($("#"+aux[0]+"_input")){
			var objId = $("#"+aux[0]+"_input").attr("data-objId");
			if (objId){
				// console.log(paramName)
				serialize_form = serialize_form.replace(paramName , aux[0]+"="+objId)
			}
		}
	});
	return serialize_form
}

function evalElements() {
	/*
	 * Function to evaluate the elements in a form, depending the type of the
	 * param, it is evaluated in a diferent way 
	 */
	
	$("tr").each(function(index) {

		// Get the identifier (id) for the parameter
		var param = $(this).attr('id');
		// Get the value for the parameter
		var value = $(this).val();
		if(value.length == 0){
			value = $(this).attr('value');
		}
		// Get the type (data-type)
		var type = $(this).attr('data-type');
		
		// DEBUG -----------------------------
		var debug_param = "PARAM:"+param;
		var debug_value = "VALUE:"+value;
		var debug_type = "TYPE:"+type;
		console.log(debug_param + "," +debug_value + "," +debug_type);
		
		// Depending of the parameter is processed
		
		if (type == "Group"){
			// Check expert level for the group content
//			var newLevel = $("select[name=expertLevel]").val();
			var newLevel = $("input[name=expertLevel]:checked").val();
			var expLevel = $(this).attr('data-expert');
			evalExpertLevel(expLevel, newLevel, $(this))
		}
		
		switch (type){
			
			case "EnumParam":
				// Get the kind of EnumParam
				var typeEnum = parseInt($(this).attr('data-enum'));
				
				if(typeEnum){
					onChangeEnumParamCombo(param + "_select", param);
				} else {
					onChangeEnumParamList(value, param);
				}
				break;
				
			case "MultiPointerParam":
				recalculateSelectDim('#' + param + '_input');
				break;
				
			default:
				onChangeParam(value, param);
				break;
		}
	});
}

function onChangeParam(value, paramId) {
	/* 
	 * Update the parameter for an element.
	 */
	setParamValue(paramId, value);
}
	

function onChangeEnumParamCombo(elemId, paramId) {
	/*
	 * Update the parameter for an element type ENUM_PARAM_COMBO.
	 */
	var elem = document.getElementById(elemId);
	setParamValue(paramId, elem.selectedIndex);
}

function onChangeEnumParamList(index, paramId) {
	/*
	 * Update the parameter for an element type ENUM_PARAM_LIST.
	 */
	setParamValue(paramId, index);
}

function setParamValue(paramId, value) {
	/*
	 * Put the new value in an attribute of the parent node.
	 */
	
	// Get the row affected
	var row = $("tr#" + paramId);
	// Update the value for the row
	row.val(value);
	
	// Get the new expertise level
//	var newLevel = $("select[name=expertLevel]").val();
	var newLevel = $("input[name=expertLevel]:checked").val();
	 
	// DEBUG
	console.log("PARAM TO EVALUATE: " + paramId)
//	console.log("WITH LEVEL: " + newLevel)
	
	// Evaluate the dependencies for the new expert level and the row affected
	evalDependencies(row, newLevel);

	// Get the params affected with the changes
	var params = row.attr('data-params');

	// Evaluate the expert level
	if (params != undefined && params.length <= 0) {
		var expLevel = row.attr('data-expert');
		evalExpertLevel(expLevel, newLevel, row)
	}
	
	
	console.log(row)
	
	// To process the hidden elements into protocol form
	// is necessary to be evaluated himself.
	evalRow(row)
}

function evalExpertLevel(expLevel, newLevel, row){
//	console.log('Evaluate the expert level')
	var expLevel = row.attr('data-expert');

	if (expLevel > newLevel) {
//		console.log("hide")
		row.hide();
	} else {
//		console.log("show")
		row.show();			
	}
}

function evalRow(row){
	var evalThis = row.attr("data-cond")
	
	console.log("EVALUATE: "+ evalThis)
	
	
	switch (evalThis){
		case "False":
			row.css('display', 'none')
			break;
		case "True":
			row.css('display', 'table-row')
			break;
	}
}


function evalDependencies(row, newLevel) {
	/*
	 * Function to evaluate the parameters dependencies for a expertise level.
	 */
	
	// Get dependencies for the parameter
	var dependencies = row.attr('data-depen');
	
//	console.log("Dependencies:", dependencies)

	if (dependencies != undefined && dependencies.length > 0) {
		
		// Dependencies splitted to be looked over the elements
		var arrayDepends = dependencies.split(",");

		for (var cont = 0; cont < arrayDepends.length; cont++) {
			// Get params affected with the dependencies
			var row2 = $("tr#" + arrayDepends[cont]);
			
//			console.log("TO EVALUATE: tr#" + arrayDepends[cont])

			// Evaluate the new parameter affected
			var res = evalCondition(row2);
			
			if (res != undefined){

				// Get the expertise level for the row affected
				var expLevel = row2.attr('data-expert');
				if (res == false || expLevel > newLevel) {
					row2.hide();
				} else if (res == true) {
					row2.show();
					
					// Evaluate the dependencies for the new row affected
					evalDependencies(row2, newLevel);
				}
			}
		}
	}	
}

function evalCondition(row) {
	/*
	 * Function to evaluate a condition given a row of the form.
	 */
	
	var res = undefined;
	var cond = row.attr('data-cond');
	
//	console.log("data-cond:"+cond)
	
	if(cond != undefined){
	
		var params = row.attr('data-params');
		
//		console.log("Params:"+ params + " length:"+params.length + "\n")
		
		if (params.length > 0){
			
			var arrayParams = params.split(",");
		
			// Get value of the element with name=itemName
			var param = null;
			var value = null;
			var cond_eval = cond;
	
			for (var cont = 0; cont < arrayParams.length; cont++) {
				param = arrayParams[cont];
				
				value = $("tr#" + param).val();
				if (!value){
					value = $("tr#" + param).attr("value");
					if (!value){
						value="''";
					}
				}
				params += "param: " + param + " value: " + value + "\n";
				cond_eval = cond_eval.replace(param, value);
			}
			
//			console.log("condition: " + cond + " \nparams:\n" + params + "\n eval: " + cond_eval);
			
			cond_eval = normalizeConditions(cond_eval)
			
			//	To check a good eval
//			console.log(cond_eval + "/" + eval(cond_eval))
		
			res = eval(cond_eval);
		}
	}
	return res;
}


function normalizeConditions(cond){
	/*
	 * For some conditions not normalized, it is replaced to be evaluated right.
	 */
	cond = replaceAll("not","!", cond);
	cond = replaceAll("and","&&", cond);
	cond = replaceAll("or","||", cond);
	cond = replaceAll("'0'","false", cond);
	cond = replaceAll("'1'","true", cond);
	return cond;
}


function browseObjects(paramName, type_param, value_param, pointerCondition, maxNumObjects) {
	/*
	 * Browse object in the database.
	 * Params: objClass: the class to get instances from (also subclasses)
	 * protClassName: class refered to a protocol
	 */
	
	 var url_param = ""
	 
     switch (type_param){
    
    	case "objClass":
			url_param = "/browse_objects/?"
				+ "&objClass=" + value_param 
				+ "&objFilter=" + pointerCondition
    		break;
    	
    	case "protClassName":
			url_param = "/browse_protocol_class/?"
				+ "&protClassName=" + value_param
			break;
			
    	case "relationClassName":
    		var res = value_param.split(",")
    		
			url_param = "/browse_relations/?"
				+ "&relationName=" + res[0]
				+ "&attributeName=" + res[1]
				+ "&protId=" + res[2]
				+ "&direction=" + res[3]
				                      
    		break;
    }
	
//	console.log("URL:", url_param)
	
	 var URL = getSubDomainURL() + url_param
	 $.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		success : function(json) {
			// specifying a dataType of json makes jQuery pre-eval the response
			// for us
			var res = getTableFormatted(paramName, json, value_param, 1);
			var selectionFunc = "processSelectionTable"
			if (maxNumObjects == 0 || maxNumObjects > 1){
				selectionFunc = "processMultipleSelectionTable"
			}
			selectDialog(paramName, res, selectionFunc);
		}
	});
}

function formProtSimple(param){
	/*
	 * Launch a custom protocol form with less options, thought for workflows
	 * where some options not need to be chosen.
	 */
	var protSimple = $("#"+param +"_input").val();
	var dataProt = $("#"+param+"_input").attr("data-prot")
	
	if (protSimple.length > 0){
		if(dataProt != undefined){
			// load the protocol params in the form
			var url = '/form/?protocolClass='+protSimple+'&action=protSimple&paramProt='+param+'&'+dataProt
		} else {
			// load a blank form with a new protocol
			var url = '/form/?protocolClass='+protSimple+'&action=protSimple&paramProt='+param
		}
		customPopup(url,500,450);
	}
	else{
		errorPopup("Error", "Protocol was not selected, please choose one.")
	}
}

function returnProtocol(){
	/*
	 * Function to return a simple protocol form with the values serializes to be
	 * saved, used for the custom workflow where the parameters need to be changed.
	 */
	params = $("#protocolForm").serialize();
	paramProt = $("#paramProt").val();
	window.opener.setParamProt(paramProt, params);
	infoPopup("Successful", "Protocol saved inside the workflow", 1);
}

function setParamProt(paramProt, params){
	/*
	 * Function to update a protocol serialize param inside the workflow.
	 * The serializes values are storaged inside the input respectively.
	 */
	$("#"+paramProt+"_input").attr("data-prot", params)
}


function getTableFormatted(node, json, id, previsualize) {
	/*
	 * This function return an html table with the elements available after use, 
	 * the function browse_object.
	 */
	var res = "<table class='content' style='overflow:auto' data-node='" + node
			+ "'><tr><th>Name</th><th>Description</th><th></th></tr>";
	
	var func = "";
	var first = "<a href='javascript:";
	var second = "'><i class='fa fa-eye'></i></a>";
	
	var x = 0;
	$.each(json, function(key, value) {
		// key is the param ObjId for the object
		// value is the name of the object
		if (value["type"] == "obj"){
			if(previsualize)
				var func = first + 'launchViewer("'+ key +'")' + second;
			res += "<tr id='"+ x + "' class='" + key + "' value='"
			+ value["nameId"]  + "' onclick=javascript:selTableMessi($(this)); ><td><i class='fa fa-file'></i>&nbsp;&nbsp;" 
			+ value["objParentName"] + "</td><td>"  + value["info"]+"</td><td>"+ func +"</td></tr>";
			
		} else if (value["type"] == "set"){
			if(previsualize)
				var func = first + 'launchViewer("'+ key +'")' + second;
			
			res += "<tr onclick='showHideChildsRow("+ key +")' class='" + key + "' value='" + value["nameId"] +"'>"
			res += "<td><i id='folderIco-"+ key +"' class='fa fa-folder'></i>&nbsp;&nbsp;"+value["objParentName"]+"</td><td>"+ value["info"]+"</td>";
			res += "<td>"+ func +"</td>"
			res += "</tr>"

			for(i=0; i<value['objects'].length; i++){
				var nameIdObj = value['objects'][i]["nameId"]
				var objIdObj = value['objects'][i]["objId"]
				var infoObj = value['objects'][i]["info"]
				
				var idText = value["nameId"] + " [item " + objIdObj +"]"
				var idElm = key + "::" + objIdObj

//				if(previsualize)
//					var func = first + 'launchViewer("'+ idElm +'")' + second;
				
				res += "<tr style='display:none;' data-row='"+ key +"' id='"+ x + "' class='" + idElm + "' value='"
				+ idText  + "' onclick=javascript:selTableMessi($(this));><td style='text-align:center;'><i class='fa fa-file'></i>&nbsp;&nbsp; item " 
				+ objIdObj + "</td><td>"+ infoObj +"</td>";
				
//				res += "<td>"+ func +"</td>"
				res += "</tr>"
					
				x++;
			}
		}
		x++;
	});
	
	res = res + "</table>";
	return res;
}

function showHideChildsRow(key){
	if($("i#folderIco-"+key).attr("class")== "fa fa-folder"){
		$("tr[data-row="+ key +"]").show();
		$("i#folderIco-"+key).attr("class", "fa fa-folder-open")
	}
	else {
		$("tr[data-row="+ key +"]").hide();
		$("i#folderIco-"+key).attr("class", "fa fa-folder")
	}
}


function selectDialog(objClass, msg, funcName) {
	/*
	 * This function create a popup with messi.js to choose the elements showed
	 * with the method getTableFormatted. 
	 */
	new Messi(msg, {
		title : 'Select ' + objClass,
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Select',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : funcName
		}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'fa-ban'
		} ]
	});
}

function getSelectedValue(elm){
	var res = []
	           
	var selected = $("tr#" + elm.attr('value'));
	res[0] = selected.attr('value');
	res[1] = selected.attr('class');
	
	return res;
}

function processSelectionTable(elm) {
	/*
	 * After select an element (PointerParam) showed with the function selectDialog with the 
	 * table generated by the function getTableFormatted, this selection is 
	 * storaged with this method.
	 */
	var value = getSelectedValue(elm);
	$('input#' + elm.attr('data-node') + '_input').attr('value', value[0]);
	$('input#' + elm.attr('data-node') + '_input').attr('data-objId', value[1]);
}

function processMultipleSelectionTable(elm) {
	/*
	 * After select an element (MultiPointerParam) showed with the function selectDialog with the 
	 * table generated by the function getTableFormatted, this selection is 
	 * storaged with this method.
	 */
	var value = getSelectedValue(elm);
	var exists = false;
	$('#'+elm.attr('data-node') + '_input option').each(function(){
		if (this.value == value[1]){
			exists = true;
			return false;
		}
	});
	
	if (exists){
		errorPopup("Selection Error","File already selected")
	}
	else{
		var selectElementName = '#'+ elm.attr('data-node') + '_input' 
		$(selectElementName).append('<option value='+value[1]+'>'+value[0]+'</option>')
		recalculateSelectDim(selectElementName)
	}	
	
}

function removeObjects(paramName){
	/*
	 * Remove element/s from multiple select element (it is usually associated to a MultiPointerParam)
	 */
	var selectElementName = '#'+paramName+ '_input'
	$(selectElementName + ' option:selected').each(function(index){
		$(this).remove();
		recalculateSelectDim(selectElementName)	
	});
}

function recalculateSelectDim(selectElementName){
	/*
	 * Calculate dim for select multiple element depending on the number of options.
	 * This is usually called after adding or removing an option 
	 */
	var selectElementHeight = $(selectElementName+' option').length * 20
	if (selectElementHeight>80){selectElementHeight=80}
	else if (selectElementHeight<23) {selectElementHeight=23}
	$(selectElementName).height(selectElementHeight)
}


function selTableMessi(elm) {
/*
 * Used to choose a element in the protocol form
 */
	var row = $("table.content");
	var id = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
//		rowOld.attr('style', '');
		// Fixed
		rowOld.css('background-color', '');
		rowOld.css('font-weight', '');
	}
	row.attr('value', id);
	elm.attr("selected","selected")
//	elm.attr('style', 'background-color: #F3CBCB;font-weight: bold;');
	// Fixed
	elm.css('background-color', '#F3CBCB')
	elm.css('font-weight', 'bold')
}

