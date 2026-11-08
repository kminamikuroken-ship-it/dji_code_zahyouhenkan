const form = document.querySelector('#coordinate-form');
const resultPanel = document.querySelector('#results');
const errorBox = document.querySelector('#error');

// JGD2011 / GRS80 origins of Japan's plane rectangular coordinate systems.
// Civil engineering convention: X = northing, Y = easting (metres).
const zoneDefinitions = {
  1: { lat: 33, lon: 129.5 },
  2: { lat: 33, lon: 131 },
  3: { lat: 36, lon: 132 + 10 / 60 },
  4: { lat: 33, lon: 133.5 },
  5: { lat: 36, lon: 134 + 20 / 60 },
  6: { lat: 36, lon: 136 },
  7: { lat: 36, lon: 137 + 10 / 60 },
  8: { lat: 36, lon: 138.5 },
  9: { lat: 36, lon: 139 + 50 / 60 },
  10: { lat: 40, lon: 140 + 50 / 60 },
  11: { lat: 44, lon: 140 + 15 / 60 },
  12: { lat: 44, lon: 142 + 15 / 60 },
  13: { lat: 44, lon: 144 + 15 / 60 },
  14: { lat: 26, lon: 142 },
  15: { lat: 26, lon: 127.5 },
  16: { lat: 26, lon: 124 },
  17: { lat: 26, lon: 131 },
  18: { lat: 20, lon: 136 },
  19: { lat: 26, lon: 154 },
};

const A = 6378137;
const INV_F = 298.257222101;
const F = 1 / INV_F;
const E2 = F * (2 - F);
const K0 = 0.9999;
const RAD = Math.PI / 180;

function projectionConstants() {
  const n = F / (2 - F);
  const n2 = n * n;
  const n3 = n2 * n;
  const n4 = n2 * n2;
  const n5 = n4 * n;
  const n6 = n3 * n3;
  return {
    Abar: A / (1 + n) * (1 + n2 / 4 + n4 / 64 + n6 / 256),
    alpha: [
      n / 2 - 2 * n2 / 3 + 5 * n3 / 16 + 41 * n4 / 180 - 127 * n5 / 288 + 7891 * n6 / 37800,
      13 * n2 / 48 - 3 * n3 / 5 + 557 * n4 / 1440 + 281 * n5 / 630 - 1983433 * n6 / 1935360,
      61 * n3 / 240 - 103 * n4 / 140 + 15061 * n5 / 26880 + 167603 * n6 / 181440,
      49561 * n4 / 161280 - 179 * n5 / 168 + 6601661 * n6 / 7257600,
      34729 * n5 / 80640 - 3418889 * n6 / 1995840,
      212378941 * n6 / 319334400,
    ],
    beta: [
      n / 2 - 2 * n2 / 3 + 37 * n3 / 96 - n4 / 360 - 81 * n5 / 512 + 96199 * n6 / 604800,
      n2 / 48 + n3 / 15 - 437 * n4 / 1440 + 46 * n5 / 105 - 1118711 * n6 / 3870720,
      17 * n3 / 480 - 37 * n4 / 840 - 209 * n5 / 4480 + 5569 * n6 / 90720,
      4397 * n4 / 161280 - 11 * n5 / 504 - 830251 * n6 / 7257600,
      4583 * n5 / 161280 - 108847 * n6 / 3991680,
      20648693 * n6 / 638668800,
    ],
  };
}

const PROJ = projectionConstants();

function meridionalArc(phi) {
  return PROJ.Abar * (phi - PROJ.alpha.reduce(
    (sum, coefficient, index) => sum + coefficient * Math.sin(2 * (index + 1) * phi),
    0,
  ));
}

function forwardTM(phi, lambda, phi0, lambda0) {
  const tau = Math.tan(phi);
  const tauPrime = Math.sinh(
    Math.asinh(tau)
      - Math.sqrt(E2) * Math.atanh(
        Math.sqrt(E2) * tau / Math.sqrt(1 + tau * tau),
      ),
  );
  const deltaLambda = lambda - lambda0;
  const xiPrime = Math.atan2(tauPrime, Math.cos(deltaLambda));
  const etaPrime = Math.asinh(
    Math.sin(deltaLambda) / Math.sqrt(tauPrime * tauPrime + Math.cos(deltaLambda) ** 2),
  );
  let xi = xiPrime;
  let eta = etaPrime;
  for (let j = 1; j <= PROJ.beta.length; j++) {
    xi -= PROJ.beta[j - 1] * Math.sin(2 * j * xiPrime) * Math.cosh(2 * j * etaPrime);
    eta -= PROJ.beta[j - 1] * Math.cos(2 * j * xiPrime) * Math.sinh(2 * j * etaPrime);
  }
  return {
    x: K0 * PROJ.Abar * (xi - meridionalArc(phi0) / PROJ.Abar),
    y: K0 * PROJ.Abar * eta,
  };
}

function inversePlaneRectangular(x, y, zone) {
  const origin = zoneDefinitions[zone];
  if (!origin) throw new Error('平面直角座標系を選択してください。');
  const phi0 = origin.lat * RAD;
  const lambda0 = origin.lon * RAD;
  let phi = phi0 + x / (A * K0);
  let lambda = lambda0 + y / (A * K0 * Math.cos(phi0));

  for (let i = 0; i < 12; i++) {
    const projected = forwardTM(phi, lambda, phi0, lambda0);
    const dx = x - projected.x;
    const dy = y - projected.y;
    if (Math.hypot(dx, dy) < 1e-7) break;

    const h = 1e-7;
    const pPhi = forwardTM(phi + h, lambda, phi0, lambda0);
    const pLambda = forwardTM(phi, lambda + h, phi0, lambda0);
    const a = (pPhi.x - projected.x) / h;
    const b = (pLambda.x - projected.x) / h;
    const c = (pPhi.y - projected.y) / h;
    const d = (pLambda.y - projected.y) / h;
    const determinant = a * d - b * c;
    if (!Number.isFinite(determinant) || Math.abs(determinant) < 1e-8) {
      throw new Error('入力座標を変換できません。X・Yと系番号を確認してください。');
    }
    phi += (dx * d - b * dy) / determinant;
    lambda += (a * dy - dx * c) / determinant;
  }

  const lat = phi / RAD;
  const lon = lambda / RAD;
  if (!Number.isFinite(lat) || !Number.isFinite(lon) || Math.abs(lat) > 90 || Math.abs(lon) > 180) {
    throw new Error('変換結果が有効な範囲外です。入力値と系番号を確認してください。');
  }
  return { lat, lon };
}

function findNumber(object, key) {
  if (!object || typeof object !== 'object') return null;
  if (Object.hasOwn(object, key)) {
    const value = Number(object[key]);
    if (Number.isFinite(value)) return value;
  }
  for (const value of Object.values(object)) {
    const found = findNumber(value, key);
    if (found !== null) return found;
  }
  return null;
}

async function getJSON(url) {
  const response = await fetch(url, { headers: { Accept: 'application/json' } });
  if (!response.ok) throw new Error(`HTTP ${response.status}`);
  return response.json();
}

function setText(id, value) { document.getElementById(id).textContent = value; }
function fmt(value) { return `${value.toFixed(3)} m`; }

form.addEventListener('submit', async (event) => {
  event.preventDefault();
  errorBox.hidden = true;
  const data = new FormData(form);
  const x = Number(data.get('x'));
  const y = Number(data.get('y'));
  const z = Number(data.get('z'));
  if (![x, y, z].every(Number.isFinite)) {
    errorBox.textContent = 'X・Y・Zに有効な数値を入力してください。';
    errorBox.hidden = false;
    return;
  }

  try {
    const { lat, lon } = inversePlaneRectangular(x, y, Number(data.get('zone')));
    const latText = lat.toFixed(9);
    const lonText = lon.toFixed(9);
    setText('latitude', latText);
    setText('longitude', lonText);
    setText('input-z', fmt(z));
    setText('geoid', '取得中…');
    setText('orthometric', '計算中…');
    setText('dem', '取得中…');
    resultPanel.hidden = false;
    resultPanel.scrollIntoView({ behavior: 'smooth', block: 'start' });

    const geoidUrl = new URL('https://vldb.gsi.go.jp/sokuchi/surveycalc/api02/geoid');
    geoidUrl.search = new URLSearchParams({ lat: latText, lon: lonText, outtype: 'JSON' });
    const demUrl = new URL('https://cyberjapandata2.gsi.go.jp/general/dem/scripts/getelevation.php');
    demUrl.search = new URLSearchParams({ lat: latText, lon: lonText, outtype: 'JSON' });
    const [geoidResult, demResult] = await Promise.allSettled([getJSON(geoidUrl), getJSON(demUrl)]);

    let geoid = null;
    if (geoidResult.status === 'fulfilled') geoid = findNumber(geoidResult.value, 'geoidHeight');
    if (geoid === null && geoidResult.status === 'fulfilled') geoid = findNumber(geoidResult.value, 'geoidheight');
    if (geoid !== null) {
      setText('geoid', fmt(geoid));
      setText('orthometric', fmt(z - geoid));
    } else {
      setText('geoid', '取得できませんでした');
      setText('orthometric', '—');
    }

    let dem = null;
    if (demResult.status === 'fulfilled') dem = findNumber(demResult.value, 'elevation');
    if (dem !== null && dem >= -9999) {
      const source = demResult.value.hsrc ? `（${demResult.value.hsrc}）` : '';
      setText('dem', `${fmt(dem)}${source}`);
    } else {
      setText('dem', '取得できませんでした');
    }

    const failures = [geoid === null && 'ジオイド高', dem === null && '地理院標高'].filter(Boolean);
    setText('api-note', failures.length
      ? `${failures.join('・')}を取得できませんでした。通信環境を確認するか、国土地理院のサイトでご確認ください。`
      : 'ジオイド高は国土地理院API、地理院標高はDEMから取得しました。地形データの解像度や取得状況により値が異なる場合があります。');
  } catch (error) {
    errorBox.textContent = error.message || '座標を変換できませんでした。入力値を確認してください。';
    errorBox.hidden = false;
  }
});

document.querySelectorAll('.copy').forEach((button) => {
  button.addEventListener('click', async () => {
    const value = document.getElementById(button.dataset.copy).textContent;
    try {
      await navigator.clipboard.writeText(value);
      button.textContent = 'コピー済み';
      setTimeout(() => { button.textContent = 'コピー'; }, 1200);
    } catch {
      button.textContent = '選択してコピー';
    }
  });
});
