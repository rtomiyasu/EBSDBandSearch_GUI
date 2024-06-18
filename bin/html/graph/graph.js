const id_graphInput = 'graph-input';
const id_graphRadon = 'graph-radon';
const id_graph1stDeriv = 'graph-1st-deriv';
const id_graph2ndDeriv = 'graph-2nd-deriv';
const line_bandEdge = 'band-edge';
const line_bandCenter= 'band-center';
const line_bandEdge2nd = 'band-edge-2nd';

function graph_updateInputImage(graphDivID, title, image, bandEdge, bandCenter) {
	ranges[line_bandEdge].range = [0, bandEdge.length]; // range for band_edge
	ranges[line_bandCenter].range = [bandEdge.length, bandEdge.length + bandCenter.length]; // range for band_center
	graph_update(graphDivID, title, '', '', image, bandEdge, bandCenter)
}

// EBSD、radon、1次微分、2次微分の各画像に共用
function graph_update(graphDivID, title, x_label, y_label, image, lines_edge, lines_center) {

	const isInput = (graphDivID === 'graph-input'); // EBSD画像か？
	const is2ndDeriv = (graphDivID === id_graph2ndDeriv); // 22次微分画像か？
	const isNoData = (image === noData); // 画像データ無しか？

	// 2次微分値画像ならば、lines_edgeの数の分だけrangeを設定
	if(is2ndDeriv) { ranges[line_bandEdge2nd].range = [0, lines_edge.length]; }

	//|
	//| data
	//|
	const data = [];
	const data_image = isNoData ? {} : {
		source : image.source, // 画像のソース（URLまたはデータURI）を指定 
		type: 'image', // https://plotly.com/javascript/reference/image/
		x0: isInput ? 0 : 0,
		y0: isInput ? 0 : (-image.size[1] / 2 + 0.5), // 画像の縦中心線をy原点
		dx: isInput ? 1 : 180 / image.size[0], // x成分、1画素 180deg/画像サイズ(W)
		dy: isInput ? 1 : 1, // y成分、１画素 そのまま
	};

	// バンドデータ収納場所
    let data_text = {
		x:[], // x座標
		y:[], // y座標
		text:[], // バンド番号（string）
		mode: 'text',
		textfont: { size: 12, color: '#000' },
    };
	// padding関数の定義(長さ(2+10の桁の数)の空白)
	const padding = i => '  ' + ' '.repeat(1 + Math.floor(Math.log10(i)));
	
	// EBSD画像の処理（）
	if(graphDivID === id_graphInput && lines_center !== null) {
		for(let i=0; i<lines_center.length; i++) {
			[lines_center[i].from, lines_center[i].to].forEach(p => {
				data_text.x.push(p[0]); // x座標
				data_text.y.push(p[1]); // y座標

				// paddingや<br>を加えることで、重なりや画像からのはみ出しをなくす。
				// 条件に基づいて prefix,suffix を設定（prefix バンド番号 suffix）
				// x座標 右端の場合：(prefix, suffix) = (padding, '')
				// y座標 上端の場合：(prefix, suffix) = (<br>, '')
				const prefix = p[0]>=(image.size[0]-1) ? padding(i+1) : p[1]>=(image.size[1]-1) ? '<br>' : '';
				
				// x座標　左端の場合：(prefix, suffix) = ('', padding)
				// y座標　下端の場合：(prefix, suffix) = ('', <br>)
				const suffix = p[0]<=0 ? padding(i+1) : p[1]<=0 ? '<br>' : '';
				// prefix と i+1 と suffix を結合して、data_text.text 配列に追加
				data_text.text.push(prefix + (i+1) + suffix);
			});
		}
	}

	// 2次微分画像の処理
	// エッジデータ＆バンド番号を収納
	if(graphDivID === id_graph2ndDeriv && lines_edge !== null) {
		for(let i=0; i<lines_edge.length; i++) {
			const p = lines_edge[i];
			data_text.x.push(p.from[0]);
			data_text.y.push(p.from[1]);
			data_text.text.push(padding(i+1) + (i+1));
		}
	}

	data.push(data_image);
	data.push(data_text);
	
	//|
	//| layout
	//|
	const xaxis = { // object from plotly
		title: x_label,
		titlefont: { size: 15 },
		showline: false,
		zeroline: false,
		ticks: 'outside',
//	    constrain: 'domain',  // 縦長（横長）のグラフの時に横軸（縦軸）を短くする（⇔'range'）type 'image'の場合は自動的に有効になる
	};
	const yaxis = { // object from plotly
		...xaxis, // xaxisを継承
		title: y_label,
	};
	
	// 線分描画map生成関数
	// 引数　pos: 座標データの配列。各要素は from と to を持ち、線分の始点と終点を表す。
	//       color: 線分の色
	const get_shapes = (pos, color) => pos === null ? [] : pos.map (
						p => ({ 
								x0: p.from[0], // 線分始点のx座標
								y0: p.from[1], // 線分始点のy座標
								x1: p.to[0],   // 線分終点のx座標
								y1: p.to[1],   // 線分終点のy座標
								type:'line',
								opacity:1.0,   // 透明度 (1:透明無し)
								line:{color:color, width:1} // 色, 幅
															}));
	// バンドエッジとバンドセンターのライン描画オブジェクト統合（EBSD画像のみ）
	const shapes = [...get_shapes(lines_edge, 'rgb(255,255,100)'), ...get_shapes(lines_center, 'rgb(100,255,100)')];
    const layout = {
		title: {
			text: title, // グラフのタイトルとして表示するテキスト
			font: { size: 20 }, // タイトルのフォント設定
			xref: 'paper', // タイトルの基準参照点を指定します。'paper' は全体の紙の領域を基準にすることを意味。
		},
		margin: { t:100, r:10, b:70, l:70 },
		// グラフのマージン（余白）サイズをピクセル単位で設定。t：上、r：右、b：下、l：左
		
		// グラフに画像を挿入。isNoData が真の場合にのみ画像を表示し、そうでない場合は空の配列を設定。
		images: isNoData ?
		[{
			"source": image.source,
			"xref": "x", // 画像の x 軸参照。ここでは 'x' として、x軸を基準。
			"yref": "y", // 画像の y 軸参照。ここでは 'y' として、y軸を基準
			"x": 0, // 画像の左下隅の x 座標。
			"y": 0, // 画像の左下隅の y 座標。
			"sizex": image.size[0], // 画像の幅
			"sizey": image.size[1], // 画像の高さ
			"sizing": "stretch", // 画像のサイズ設定方法。'stretch' は画像を伸縮させて指定されたサイズに合わせる。
		}] : [], // 空の配列
		hovermode: is2ndDeriv ? 'closest' : false, // falseにするとホバーが表示されなくなるが、クリックイベントが発生しなくなる
		// hovermode が 'closest' に設定されている場合、ホバーツールチップはカーソルに最も近いデータポイントの情報を表示
        xaxis, // xaix設定変数
        yaxis, // yaxis設定変数
		shapes, // グラフに描画するシェイプを指定します
    };

	// no dataの場合
	if(isNoData) {
		xaxis.range = [0, image.size[0]];
		yaxis.range = [image.size[1], 0];
	}

	// EBSD画像の場合
	else if(isInput) {
//		yaxis.scaleanchor = 'x'; // グリッドを正方形にする。type 'image'の場合は自動的に有効になる
		yaxis.constraintoward = 'top';
		// グラフの上端を基準にしてスケールが固定される。これにより、x軸とy軸のスケールが上方向に向かって調整される。
		// データポイントが上方向に寄るように表示され、下側の余白が増える。
		layout.margin = { t:100, r:10, b:20, l:40 };
	}
	else {
		yaxis.scaleratio = 180 / image.size[1];
		// radon画像の場合は、縦軸のスケールを横軸180に合わせている。
	}
	const xrange = [0, isNoData ? image.size[0] : isInput ? image.size[0]-1 : 180];
	// isNoData : [0, image.size[0]], isInput : [0, image.size[0]-1, else [0, 180]]
	const yrange = isNoData ? [image.size[1], 0] : isInput ? [image.size[1]-1, 0] : [image.size[1]/2, -image.size[1]/2];
	// isNoData : [image.size[1], 0], isInput : [image.size[1]-1, 0], else [image.size[1]/2, -image.size[1]/2] 
	const ratio = Math.min(1, (yrange[0]-yrange[1])/(xrange[1]-xrange[0])); //1を超えないようにして、縦長にならないようにしている。
	
	// 指定されたID (graphDivID) を持つHTML要素を取得
	const div = document.getElementById(graphDivID);
	// HTML要素の高さを計算
	// ビューポート幅（vw）の100倍に ratio を掛けた値に、レイアウトの上下と左右のマージンの差を調整した値を足したもの。
	// vw(Viewport Width) : 現在表示されている領域の幅
	div.style.height = `calc(${100 * ratio}vw + ${(layout.margin.t + layout.margin.b) - (layout.margin.l + layout.margin.r) * ratio}px)`;

	//|
	//| config
	//|
    const config = {
        displayModeBar: false, // 拡大縮小ボタンなどを非表示にする
        responsive: true, // ウィンドウのリサイズ時に、ウィンドウサイズに合わせて再描画する
//        staticPlot: true, // ツールチップ表示やズームを無効にする
    };

	//Plotly.newPlot(graphDiv, data, layout, config);
	//graphDiv: グラフを描画するDOM要素のIDや要素そのもの。
	//data: 描画するデータを含む配列。
	//layout: グラフのレイアウトやスタイルを定義するオブジェクト。
	//config: オプションの設定を定義するオブジェクト（省略可能）

    Plotly.newPlot(graphDivID, data, layout, config);
}


const ranges = {};
ranges[line_bandEdge] = {id: id_graphInput, range:[0,0]};
ranges[line_bandCenter] = {id: id_graphInput, range:[0,0]};
ranges[line_bandEdge2nd] = {id: id_graph2ndDeriv, range:[0,0]};
const noData = { size:[100,80], source:'data:image/svg+xml;charset=utf8,%3C%3Fxml%20version%3D%221.0%22%20encoding%3D%22UTF-8%22%3F%3E%0A%3Csvg%20width%3D%22100mm%22%20height%3D%2280mm%22%20version%3D%221.1%22%20viewBox%3D%220%200%20100%2080%22%20xmlns%3D%22http%3A%2F%2Fwww.w3.org%2F2000%2Fsvg%22%20xmlns%3Acc%3D%22http%3A%2F%2Fcreativecommons.org%2Fns%23%22%20xmlns%3Adc%3D%22http%3A%2F%2Fpurl.org%2Fdc%2Felements%2F1.1%2F%22%20xmlns%3Ardf%3D%22http%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23%22%3E%3Cmetadata%3E%3Crdf%3ARDF%3E%3Ccc%3AWork%20rdf%3Aabout%3D%22%22%3E%3Cdc%3Aformat%3Eimage%2Fsvg%2Bxml%3C%2Fdc%3Aformat%3E%3Cdc%3Atype%20rdf%3Aresource%3D%22http%3A%2F%2Fpurl.org%2Fdc%2Fdcmitype%2FStillImage%22%2F%3E%3Cdc%3Atitle%2F%3E%3C%2Fcc%3AWork%3E%3C%2Frdf%3ARDF%3E%3C%2Fmetadata%3E%3Cg%20transform%3D%22translate(-.25935%20-7.5665)%22%20fill%3D%22%230019a4%22%3E%3Cpath%20d%3D%22m20.509%2054.109h-2.1533l-6.205-11.707v11.707h-1.626v-13.087h2.6982l5.6601%2010.687v-10.687h1.626z%22%2F%3E%3Cpath%20d%3D%22m33.745%2042.526q0.7998%200.8789%201.2217%202.1533%200.43066%201.2744%200.43066%202.8916t-0.43945%202.9004q-0.43066%201.2744-1.2129%202.1269-0.80858%200.88768-1.916%201.3359-1.0986%200.44824-2.5136%200.44824-1.3799%200-2.5136-0.45703-1.125-0.45703-1.916-1.3271-0.79101-0.87011-1.2217-2.1357-0.42187-1.2656-0.42187-2.8916%200-1.5996%200.42187-2.8652%200.42187-1.2744%201.2305-2.1797%200.77343-0.86132%201.916-1.3183%201.1514-0.45703%202.5049-0.45703%201.4062%200%202.5224%200.46582%201.125%200.45703%201.9072%201.3096zm-0.1582%205.0449q0-2.5488-1.1426-3.9287-1.1426-1.3887-3.1201-1.3887-1.9951%200-3.1377%201.3887-1.1338%201.3799-1.1338%203.9287%200%202.5752%201.1601%203.9462%201.1601%201.3623%203.1113%201.3623t3.1025-1.3623q1.1601-1.3711%201.1601-3.9462z%22%2F%3E%3Cpath%20d%3D%22m55.542%2047.579q0%201.7842-0.78222%203.2343-0.77343%201.4502-2.0654%202.25-0.89648%200.5537-2.0039%200.7998-1.0986%200.24609-2.9004%200.24609h-3.3046v-13.087h3.2695q1.916%200%203.041%200.28125%201.1338%200.27246%201.916%200.75585%201.3359%200.83495%202.083%202.2236%200.74706%201.3887%200.74706%203.2959zm-1.8193-0.02637q0-1.5381-0.53613-2.5927-0.53613-1.0547-1.5996-1.6611-0.77343-0.43945-1.6435-0.60644-0.87011-0.17578-2.083-0.17578h-1.6347v10.099h1.6347q1.2568%200%202.1885-0.18457%200.94042-0.18457%201.7226-0.68554%200.97558-0.62402%201.459-1.6435%200.49218-1.0195%200.49218-2.5488z%22%2F%3E%3Cpath%20d%3D%22m68.69%2054.109h-1.8545l-1.2832-3.6474h-5.6601l-1.2832%203.6474h-1.7666l4.7636-13.087h2.3203zm-3.6738-5.1415-2.2939-6.4247-2.3027%206.4247z%22%2F%3E%3Cpath%20d%3D%22m78.955%2042.57h-4.6757v11.54h-1.7402v-11.54h-4.6757v-1.5469h11.092z%22%2F%3E%3Cpath%20d%3D%22m89.994%2054.109h-1.8545l-1.2832-3.6474h-5.6601l-1.2832%203.6474h-1.7666l4.7636-13.087h2.3203zm-3.6738-5.1415-2.2939-6.4247-2.3027%206.4247z%22%2F%3E%3C%2Fg%3E%3C%2Fsvg%3E%0A' };

// Plotly.restyleの入力変数：Plotly.restyle (graphDiv, update, traces)
// graphDiv : 描画されているHTML要素、update：更新する内容、trace：更新の対象となるトレース番号
function checkChanged(elem) {
	// ラベル数字の色
	if (elem.id === 'band-number-is-black-2nd') {
		Plotly.restyle(id_graph2ndDeriv, { 'textfont.color': elem.checked ? '#000' : '#FFF' }, 1); // #000：黒、#FFF：白
	}
	// ラベル数字の表示・非表示
	else if (elem.id.search('number') !== -1) {
		// elem.cheked ならば、visibleにする。
		if (elem.id === 'band-number') { Plotly.restyle(id_graphInput, { visible: elem.checked }, 1); } // 1 : 更新するトレース番号
		if (elem.id === 'band-number-2nd') { Plotly.restyle(id_graph2ndDeriv, { visible: elem.checked }, 1); }
	}
	// バンド境界・中心の表示・非表示 (EBSD画像：'band','center', 2次微分画像: 'band')
	else {
		const graphDivID = ranges[elem.id].id; // id_graphInput or id_graph2ndDeriv
		// elem_id : line_bandEdge, line_bandCenter, line_bandEdge2nd
		// id : id_graphInput, id_graph2ndDeriv
		const range = ranges[elem.id].range; // band indices range [min, max]
		// elem.idが、id_graphInputの時には、
		// 'band', 'center'のチェックに応じて、shapeが選択され、表示・非表示が切り替わる。
		const target = {};
		for(let i=range[0]; i<range[1]; i++) {
			target['shapes[' + i + '].visible'] = elem.checked;
		} // 描画するshapes(線分)を選択
		Plotly.relayout(graphDivID, target);
	}
}

function graph_updateAll(data, lines) {
    graph_updateInputImage(id_graphInput, 'Input image',
						data[id_graphInput], lines[line_bandEdge], lines[line_bandCenter]);
    graph_update(id_graph2ndDeriv, 'Second derivative', 'Projection angle (deg)', 'Projection position (pixels)',
						data[id_graph2ndDeriv], lines[line_bandEdge2nd], null);
    graph_update(id_graph1stDeriv, 'First derivative', 'Projection angle (deg)', 'Projection position (pixels)',
						data[id_graph1stDeriv], null, null);
    graph_update(id_graphRadon, 'Radon tranform (sinogram)', 'Projection angle (deg)', 'Projection position (pixels)',
						data[id_graphRadon], null, null);
	Array.from(document.querySelectorAll('input[type=checkbox]')).map(e => checkChanged(e));
	// Array.fromは、配列ライクなオブジェクト（例えばNodeList）やイテラブルオブジェクトを本物の配列に変換するメソッド。
    // ここでは、NodeListを配列に変換している。
	// document.querySelectorAllは、CSSセレクターを使って指定された条件に一致するすべての要素を取得する。
	// この場合、'input[type=checkbox]'というセレクターを使って、ページ内のすべてのチェックボックス
	// （<input>要素でtype属性がcheckboxのもの）を選択している。querySelectorAllは、NodeList（ノードリスト）を返す。
	// mapメソッドは、配列内のすべての要素に対して指定された関数を呼び出し、その結果を新しい配列として返す。
	// ここでは、配列内の各チェックボックス要素（e）に対して、checkChanged(e)関数を呼び出している。
	// checkChangedは、チェックボックスの要素eに対して、それぞれの処理を行う。(band_center表示or非表示など)

	const graph2nd = document.getElementById(id_graph2ndDeriv); // 2次微分画像の要素を取得
	// plotlyグラフのクリックイベント (plotly_click) に対するイベントハンドラを設定
	// このイベントハンドラは、グラフがクリックされたときに呼び出され、クリックに関する情報が data オブジェクトとして渡される。
    graph2nd.on('plotly_click', function(data){ 
		if(data.event.button === 2) { // クリックされたマウスボタン,2 : 右クリック
			notify(`add ${data.points[0].x} ${data.points[0].y}`);
			// data.points[0].x と data.points[0].y は、クリックされたデータ点のx座標とy座標
		}
	});
    graph2nd.on('plotly_doubleclick', data => Plotly.Fx.hover(id_graph2ndDeriv, []));
	graph2nd.on('plotly_hover', data => Plotly.Fx.hover(id_graph2ndDeriv, [])); // ホバーを非表示にする（ただしクリックすると表示される）
}

function graph_init() {
	const data = {};
	data[id_graphInput] = noData;
	data[id_graph2ndDeriv] = noData;
	data[id_graph1stDeriv] = noData;
	data[id_graphRadon] = noData;
	const lines = {};
	lines[line_bandEdge] = [];
	lines[line_bandCenter] = [];
	lines[line_bandEdge2nd] = [];
	graph_updateAll(data, lines);
}