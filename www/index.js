import * as wasm from "./dicodon-optimizer-wasm";

wasm.start();
//document.getElementById("output").innerHTML = wasm.get_back('shu');
document.getElementById("loading").style.display = 'none';
document.getElementById("interface").style.display = 'block';
//document.getElementById("dicodon_frequencies").value = '';


function btn_count(ev) {
	var file = document.getElementById('input_file');
	if (file.files.length == 0){
		alert("Pick a file first");
		return;
	}

	document.getElementById("dicodon_frequencies").value = 'calculating...';
    if(file.files.length)
    {
        var reader = new FileReader();

        reader.onload = function(e)
        {
			var res = wasm.fasta_to_score_table(e.target.result);
			res = "# from " + file.files[0].name + "\n" + res;
            document.getElementById('dicodon_frequencies').value = res;
        };

        reader.readAsBinaryString(file.files[0]);
    }
}

function btn_analyze(ev) {
    var dicodon_freqs = document.getElementById("dicodon_frequencies").value;
	if (dicodon_freqs.length == 0) {
		alert("create/enter a dicodon frequency table first");
	}
	dicodon_freqs = dicodon_freqs.replace(/[\t ]+/g,"\t");
	var seq = document.getElementById("input_sequence").value;
	if (seq.length == 0) {
		alert("enter a sequence first to analyze");
	}
	document.getElementById('output').innerHTML = "calculating";
	try {
		var output = wasm.analyze(dicodon_freqs, seq);
		var htmlout = "<table class='result_table'><tr><th>p</th><th>Codon</th><th>AA</th><th>p</th></td>";
		for (var ii in output.codons) {
			htmlout += "<tr>";
			var pout = "";
			var side = "";
			if (ii % 2 == 0) {
				side = "left";
			}
			else {
				side="right";
			}

			if (ii < output.scores.length) {
				let score = parseFloat(output.scores[ii]);
				var score_class = '';
				if (score >= 0.5) {
					score_class = "marker_ok";

				} else if (score >= 0.1) {
					score_class = "marker_passable";
				}
				else {
					score_class = 'marker_low';
				}

				pout = "<td rowspan='2' style='vertical-align:center' class='" + score_class + " " + side +"'>" + score.toFixed(2) + "</td>"
			} else {
				pout = "<td></td>";
			};

			if (ii % 2 == 0) {
				htmlout += pout;
			}else {
				htmlout += "";
			};

			htmlout += "<td>" + output.codons[ii] + "</td>"
			htmlout += "<td>" + output.amino_acids[ii] + "</td>"
			if (ii % 2 != 0) {
				htmlout += pout;
			} else {
				htmlout += "";
			};
			htmlout += "</tr>";
		}
		htmlout += "</table>"
		document.getElementById('output').innerHTML = htmlout;
		document.getElementById('output_sequence').value = output.codons.join('');
	} catch(error) {
		document.getElementById('output').innerHTML = "<span class='error'>" + error + "</span>";
		document.getElementById('output_sequence').value = '';
	}

	document.getElementById('output').scrollIntoView();
}
function btn_optimize(ev) {
	alert("todo");
}

document.getElementById("btn_count").onclick = btn_count;
document.getElementById("btn_analyze").onclick = btn_analyze;
document.getElementById("btn_optimize").onclick = btn_optimize;
