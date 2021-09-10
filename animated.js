var Windy = function( params ){
  var VELOCITY_SCALE = 0.25;        //速度比例
  var INTENSITY_SCALE_STEP = 0.1;  // 粒子强度色标的步长
  var MAX_WIND_INTENSITY = 2;      //粒子强度最大时的速度
  var MAX_PARTICLE_AGE = 10;       //在再生之前绘制粒子的最大帧数
  var PARTICLE_LINE_WIDTH = 0.5;   //绘制粒子的线宽
  var PARTICLE_MULTIPLIER = 6/100  //粒子计数标量
  var PARTICLE_REDUCTION = 0.1;    //将粒子数减少到移动设备的正常水平 
  var FRAME_RATE = 60;             //每帧所需的毫秒数 
  var NULL_WIND_VECTOR = [NaN, NaN, null];
  // 双线性内插
  var bilinearInterpolateVector = function(x, y, g00, g10, g01, g11) {
      var rx = (1 - x);
      var ry = (1 - y);
      var a = rx * ry,  b = x * ry,  c = rx * y,  d = x * y;
      var u = g00[0] * a + g10[0] * b + g01[0] * c + g11[0] * d;
      var v = g00[1] * a + g10[1] * b + g01[1] * c + g11[1] * d;
      return [u, v, Math.sqrt(u * u + v * v)];
  };
  var createWindBuilder = function(uComp, vComp) {
      var uData = uComp.data, vData = vComp.data;
      return {
          header: uComp.header,
          data: function(i) {
              return [uData[i], vData[i]];
          },
          interpolate: bilinearInterpolateVector
      }
  };
  var createBuilder = function(data) {
      var uComp = null, vComp = null, scalar = null;
      data.forEach(function(record) {
          switch (record.header.parameterCategory + "," + record.header.parameterNumber) {
              case "2,2": uComp = record; break;
              case "2,3": vComp = record; break;
              default:
                scalar = record;
          }
      });
      return createWindBuilder(uComp, vComp);
  };
  var buildGrid = function(data, callback) {
      var builder = createBuilder(data);
      var header = builder.header;
      var λ0 = header.lo1, φ0 = header.la1;  //遍历原点
      var Δλ = header.dx, Δφ = header.dy;    //遍历间隔
      var ni = header.nx, nj = header.ny;    //每行每列遍历数量
      var date = new Date(header.refTime);
      date.setHours(date.getHours() + header.forecastTime);

      //假设从0开始0。经度从 λ0 开始增加，纬度从 φ0 开始减少。
      var grid = [], p = 0;
      var isContinuous = Math.floor(ni * Δλ) >= 360;
      for (var j = 0; j < nj; j++) {
          var row = [];
          for (var i = 0; i < ni; i++, p++) {
              row[i] = builder.data(p);
          }
          if (isContinuous) {
              //复制第一列作为最后一列以简化插值逻辑
              row.push(row[0]);
          }
          grid[j] = row;
      }
      function interpolate(λ, φ) {
          var i = floorMod(λ - λ0, 360) / Δλ;
          var j = (φ0 - φ) / Δφ;
          var fi = Math.floor(i), ci = fi + 1;
          var fj = Math.floor(j), cj = fj + 1;
          var row;
          if ((row = grid[fj])) {
              var g00 = row[fi];
              var g10 = row[ci];
              if (isValue(g00) && isValue(g10) && (row = grid[cj])) {
                  var g01 = row[fi];
                  var g11 = row[ci];
                  if (isValue(g01) && isValue(g11)) {
                      // 找到所有四个点进行插值。
                      return builder.interpolate(i - fi, j - fj, g00, g10, g01, g11);
                  }
              }
          }
          return null;
      }
      callback( {
          date: date,
          interpolate: interpolate
      });
  };
  var isValue = function(x) {
      return x !== null && x !== undefined;
  }
  var floorMod = function(a, n) {
      return a - n * Math.floor(a / n);
  }
  var clamp = function(x, range) {
      return Math.max(range[0], Math.min(x, range[1]));
  }
  var isMobile = function() {
      return (/android|blackberry|iemobile|ipad|iphone|ipod|opera mini|webos/i).test(navigator.userAgent);
  }

  //计算 (x, y) 点投影形状引起的风矢量失真
  var distort = function(projection, λ, φ, x, y, scale, wind, windy) {
      var u = wind[0] * scale;
      var v = wind[1] * scale;
      var d = distortion(projection, λ, φ, x, y, windy);
      // 通过 u 和 v 缩放失真向量，然后相加。
      wind[0] = d[0] * u + d[2] * v;
      wind[1] = d[1] * u + d[3] * v;
      return wind;
  };
  var distortion = function(projection, λ, φ, x, y, windy) {
      var τ = 2 * Math.PI;
      var H = Math.pow(10, -5.2);
      var hλ = λ < 0 ? H : -H;
      var hφ = φ < 0 ? H : -H;

      var pλ = project(φ, λ + hλ,windy);
      var pφ = project(φ + hφ, λ, windy);
      //子午线比例因子其中 R = 1。这处理长度为 1º λ 的问题根据 φ 变化。否则，两极就会产生挤压效应。 
      var k = Math.cos(φ / 360 * τ);
      return [
          (pλ[0] - x) / hλ / k,
          (pλ[1] - y) / hλ / k,
          (pφ[0] - x) / hφ,
          (pφ[1] - y) / hφ
      ];
  };
  var createField = function(columns, bounds, callback) {
      function field(x, y) {
          var column = columns[Math.round(x)];
          return column && column[Math.round(y)] || NULL_WIND_VECTOR;
      }
      field.release = function() {
          columns = [];
      };
      field.randomize = function(o) {
          var x, y;
          var safetyNet = 0;
          do {
              x = Math.round(Math.floor(Math.random() * bounds.width) + bounds.x);
              y = Math.round(Math.floor(Math.random() * bounds.height) + bounds.y)
          } while (field(x, y)[2] === null && safetyNet++ < 30);
          o.x = x;
          o.y = y;
          return o;
      };
      callback( bounds, field );
  };
  var buildBounds = function( bounds, width, height ) {
      var upperLeft = bounds[0];
      var lowerRight = bounds[1];
      var x = Math.round(upperLeft[0]); 
      var y = Math.max(Math.floor(upperLeft[1], 0), 0);
      var xMax = Math.min(Math.ceil(lowerRight[0], width), width - 1);
      var yMax = Math.min(Math.ceil(lowerRight[1], height), height - 1);
      return {x: x, y: y, xMax: width, yMax: yMax, width: width, height: height};
  };
  var deg2rad = function( deg ){
    return (deg / 180) * Math.PI;
  };
  var rad2deg = function( ang ){
    return ang / (Math.PI/180.0);
  };
  var invert = function(x, y, windy){
    var mapLonDelta = windy.east - windy.west;
    var worldMapRadius = windy.width / rad2deg(mapLonDelta) * 360/(2 * Math.PI);
    var mapOffsetY = ( worldMapRadius / 2 * Math.log( (1 + Math.sin(windy.south) ) / (1 - Math.sin(windy.south))  ));
    var equatorY = windy.height + mapOffsetY;
    var a = (equatorY-y)/worldMapRadius;
    var lat = 180/Math.PI * (2 * Math.atan(Math.exp(a)) - Math.PI/2);
    var lon = rad2deg(windy.west) + x / windy.width * rad2deg(mapLonDelta);
    return [lon, lat];
  };
  var mercY = function( lat ) {
    return Math.log( Math.tan( lat / 2 + Math.PI / 4 ) );
  };
  var project = function( lat, lon, windy) { //均以弧度表示
    var ymin = mercY(windy.south);
    var ymax = mercY(windy.north);
    var xFactor = windy.width / ( windy.east - windy.west );
    var yFactor = windy.height / ( ymax - ymin );
    var y = mercY( deg2rad(lat) );
    var x = (deg2rad(lon) - windy.west) * xFactor;
    var y = (ymax - y) * yFactor; // y 指向南 
    return [x, y];
  };
  var interpolateField = function( grid, bounds, extent, callback ) {
    var projection = {};
    var velocityScale = VELOCITY_SCALE;
    var columns = [];
    var x = bounds.x;
    function interpolateColumn(x) {
        var column = [];
        for (var y = bounds.y; y <= bounds.yMax; y += 2) {
                var coord = invert( x, y, extent );
                if (coord) {
                    var λ = coord[0], φ = coord[1];
                    if (isFinite(λ)) {
                        var wind = grid.interpolate(λ, φ);
                        if (wind) {
                            wind = distort(projection, λ, φ, x, y, velocityScale, wind, extent);
                            column[y+1] = column[y] = wind;
                        }
                    }
                }
        }
        columns[x+1] = columns[x] = column;
    }
    (function batchInterpolate() {
                var start = Date.now();
                while (x < bounds.width) {
                    interpolateColumn(x);
                    x += 2;
                    if ((Date.now() - start) > 1000) {
                        setTimeout(batchInterpolate, 25);
                        return;
                    }
                }
          createField(columns, bounds, callback);
    })();
  };
  var animate = function(bounds, field) {
    function asColorStyle(r, g, b, a) {
        return "rgba(" + 243 + ", " + 243 + ", " + 238 + ", " + a + ")";
    }
    function hexToR(h) {return parseInt((cutHex(h)).substring(0,2),16)}
    function hexToG(h) {return parseInt((cutHex(h)).substring(2,4),16)}
    function hexToB(h) {return parseInt((cutHex(h)).substring(4,6),16)}
    function cutHex(h) {return (h.charAt(0)=="#") ? h.substring(1,7):h}

    function windIntensityColorScale(step, maxWind) {

        var result = [
          "rgba(" + hexToR('#006EFF') + ", " + hexToG('#006EFF') + ", " + hexToB('#006EFF') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#0084FF') + ", " + hexToG('#0084FF') + ", " + hexToB('#0084FF') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#00A2FF') + ", " + hexToG('#00A2FF') + ", " + hexToB('#00A2FF') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#00C3FFF5') + ", " + hexToG('#00C3FFF5') + ", " + hexToB('#00C3FFF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#3CF8FFF5') + ", " + hexToG('#3CF8FFF5') + ", " + hexToB('#3CF8FFF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#64FFEAF5') + ", " + hexToG('#64FFEAF5') + ", " + hexToB('#64FFEAF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#64FFC4F5') + ", " + hexToG('#64FFC4F5') + ", " + hexToB('#64FFC4F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#88FFCDF5') + ", " + hexToG('#88FFCDF5') + ", " + hexToB('#88FFCDF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#C7FF86F5') + ", " + hexToG('#C7FF86F5') + ", " + hexToB('#C7FF86F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#DCFF8AF5') + ", " + hexToG('#DCFF8AF5') + ", " + hexToB('#DCFF8AF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#E3FF43F5') + ", " + hexToG('#E3FF43F5') + ", " + hexToB('#E3FF43F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FFF461F5') + ", " + hexToG('#FFF461F5') + ", " + hexToB('#FFF461F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FFE053F5') + ", " + hexToG('#FFE053F5') + ", " + hexToB('#FFE053F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FFB350F5') + ", " + hexToG('#FFB350F5') + ", " + hexToB('#FFB350F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FF8944F5') + ", " + hexToG('#FF8944F5') + ", " + hexToB('#FF8944F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FF7D4AF5') + ", " + hexToG('#FF7D4AF5') + ", " + hexToB('#FF7D4AF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FC6449F5') + ", " + hexToG('#FC6449F5') + ", " + hexToB('#FC6449F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FF5F53F5') + ", " + hexToG('#FF5F53F5') + ", " + hexToB('#FF5F53F5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FD4E4EF5') + ", " + hexToG('#FD4E4EF5') + ", " + hexToB('#FD4E4EF5') + ", " + 0.25 + ")",
          "rgba(" + hexToR('#FF5555F5') + ", " + hexToG('#FF5555F5') + ", " + hexToB('#FF5555F5') + ", " + 0.25 + ")"
        ]
        result.indexFor = function(m) {  //将速度映射到样式 
            return Math.floor(Math.min(m, maxWind) / maxWind * (result.length - 1));
        };
        return result;
    }
    var colorStyles = windIntensityColorScale(INTENSITY_SCALE_STEP, MAX_WIND_INTENSITY);
    var buckets = colorStyles.map(function() { return []; });
    var particleCount = Math.round(bounds.width * bounds.height * PARTICLE_MULTIPLIER);
    if (isMobile()) {
      particleCount *= PARTICLE_REDUCTION;
    }
    var fadeFillStyle = "rgba(0, 0, 0, 0.97)";
    var particles = [];
    for (var i = 0; i < particleCount; i++) {
        particles.push(field.randomize({age: Math.floor(Math.random() * MAX_PARTICLE_AGE) + 0}));
    }
    function evolve() {
        buckets.forEach(function(bucket) { bucket.length = 0; });
        particles.forEach(function(particle) {
            if (particle.age > MAX_PARTICLE_AGE) {
                field.randomize(particle).age = 0;
            }
            var x = particle.x;
            var y = particle.y;
            var v = field(x, y);  //当前位置的向量 
            var m = v[2];
            if (m === null) {
                particle.age = MAX_PARTICLE_AGE;  // 粒子已经逃离网格
            }
            else {
                var xt = x + v[0];
                var yt = y + v[1];
                if (field(xt, yt)[2] !== null) {
                    // 从 (x,y) 到 (xt,yt) 的路径是可见的
                    particle.xt = xt;
                    particle.yt = yt;
                    buckets[colorStyles.indexFor(m)].push(particle);
                }
                 else {
                    particle.x = xt;
                    particle.y = yt;
                }
            }
            particle.age += 1;
        });
    }
    var g = params.canvas.getContext("2d");
    g.lineWidth = PARTICLE_LINE_WIDTH;
    g.fillStyle = fadeFillStyle;
    function draw() {
        //淡化现有的粒子轨迹
        var prev = g.globalCompositeOperation;
        g.globalCompositeOperation = "destination-in";
        g.fillRect(bounds.x, bounds.y, bounds.width, bounds.height);
        g.globalCompositeOperation = prev;
        //绘制新的粒子轨迹
        buckets.forEach(function(bucket, i) {
            if (bucket.length > 0) {
                g.beginPath();
                g.strokeStyle = colorStyles[i];
                bucket.forEach(function(particle) {
                    g.moveTo(particle.x, particle.y);
                    g.lineTo(particle.xt, particle.yt);
                    particle.x = particle.xt;
                    particle.y = particle.yt;
                });
                g.stroke();
            }
        });
    }
    (function frame() {
        try {
            windy.timer = setTimeout(function() {
              requestAnimationFrame(frame);
              evolve();
              draw();
            }, 1000 / FRAME_RATE);
        }
        catch (e) {
            console.error(e);
        }
    })();
  }
  var start = function( bounds, width, height, extent ){
    var mapBounds = {
      south: deg2rad(extent[0][1]),
      north: deg2rad(extent[1][1]),
      east: deg2rad(extent[1][0]),
      west: deg2rad(extent[0][0]),
      width: width,
      height: height
    };
    stop();
    // 构建网格 
    buildGrid( params.data, function(grid){
      //插值字段 
      interpolateField( grid, buildBounds( bounds, width, height), mapBounds, function( bounds, field ){
        //用随机点画布 
        windy.field = field;
        animate( bounds, field );
      });
    });
  };
  var stop = function(){
    if (windy.field) windy.field.release();
    if (windy.timer) clearTimeout(windy.timer)
  };
  var windy = {
    params: params,
    start: start,
    stop: stop
  };
  return windy;
}
window.requestAnimationFrame = (function(){
  return  window.requestAnimationFrame       ||
          window.webkitRequestAnimationFrame ||
          window.mozRequestAnimationFrame    ||
          window.oRequestAnimationFrame ||
          window.msRequestAnimationFrame ||
          function( callback ){
            window.setTimeout(callback, 1000 / 20);
          };
})();
