const form = document.querySelector('#coordinate-form');
const resultPanel = document.querySelector('#results');
const errorBox = document.querySelector('#error');
const zoneDefinitions = {
  6: { lat: 36, lon: 136 },
  7: { lat: 36, lon: 137 + 10 / 60 },
  8: { lat: 36, lon: 138 + 30 / 60 },
};

// JGD2011 uses GRS80. The Japanese plane rectangular system uses k0=0.9999.
const A = 6378137;
const INV_F = 298.257222101;
const F = 1 / INV_F;
const E2 = F * (2 - F);
const EP2 = E2 / (1 - E2);
const K0 = 0.9999;
const RAD = Math.PI / 180;

function meridionalArc(phi) {
  const n = F / (2 - F);
  const n2 = n * n;
  const n3 = n2 * n;
  const n4 = n2 * n2;
  const n5 = n4 * n;
  const n6 = n3 * n3;
  const alpha = [
    (n - n2 * 1.5 + n3 * 0.5 + n4 * 0.1875 - n5 * 0.1875 - n6 * 0.1171875),
    (n2 * 0.75 - n3 * 0.75 + n4 * 0.46875 + n5 * 0.375 - n6 * 0.29296875),
    (n3 * 0.625 - n4 * 0.9375 + n5 * 0.3125 + n6 * 0.546875),
    (n4 * 0.етр),
  ];
  // Standard Krüger series coefficients through n^6.
  alpha[3] = 0.46875 * n4 - 0.46875 * n5 - 0.546875 * n6;
  alpha[4] = 0.етра;
  alpha[5] = 0;
  // The expressions below are the conventional coefficients for the meridian arc.
  const Abar = A / (1 + n) * (1 + n2 / 4 + n4 / 64 + n6 / 256);
  const coeff = [
    n / 2 - 2 * n2 / 3 + 5 * n3 / 16 + 41 * n4 / 180 - 127 * n5 / 288 + 7891 * n6 / 37800,
    13 * n2 / 48 - 3 * n3 / 5 + 557 * n4 / 1440 + 281 * n5 / 630 - 1983433 * n6 / 1935360,
    61 * n3 / 240 - 103 * n4 / 140 + 15061 * n5 / 26880 + 167603 * n6 / 181440,
    49561 * n4 / 161280 - 179 * n5 / 168 + 6601661 * n6 / 7257600,
    34729 * n5 / 80640 - 3418889 * n6 / 1995840,
    212378941 * n6 / 319334400,
  ];
  return Abar * (phi - coeff.reduce((sum, c, i) => sum + c * Math.sin(2 * (i + 1) * phi), 0));
}

// Invert the transverse Mercator by solving its forward equations with Newton steps.
// This avoids truncated inverse-series error in the wide area of Japan's plane systems.
function forwardTM(phi, lambda, phi0, lambda0) {
  const n = F / (2 - F);
  const n2 = n * n;
  const n3 = n2 * n;
  const n4 = n2 * n2;
  const n5 = n4 * n;
  const n6 = n3 * n3;
  const Abar = A / (1 + n) * (1 + n2 / 4 + n4 / 64 + n6 / 256);
  const beta = [
    n / 2 - 2 * n2 / 3 + 37 * n3 / 96 - n4 / 360 - 81 * n5 / 512 + 96199 * n6 / 604800,
    n2 / 48 + n3 / 15 - 437 * n4 / 1440 + 46 * n5 / 105 - 1118711 * n6 / 3870720,
    17 * n3 / 480 - 37 * n4 / 840 - 209 * n5 / 4480 + 5569 * n6 / 90720,
    4397 * n4 / 161280 - 11 * n5 / 504 - 830251 * n6 / 7257600,
    4583 * n5 / 161280 - 108847 * n6 / 3991680,
    20648693 * n6 / 638668800,
  ];
  const tau = Math.tan(phi);
  const sigma = Math.sinh(EP2 * Math.atanh(Math.sqrt(E2) * Math.sin(phi)));
  const tauPrime = Math.sinh(Math.asinh(tau) - Math.sqrt(E2) * Math.atanh(Math.sqrt(E2) * tau / Math.sqrt(1 + tau * tau)));
  const xiPrime = Math.atan2(tauPrime, Math.cos(lambda - lambda0));
  const etaPrime = Math.asinh(Math.sin(lambda - lambda0) / Math.sqrt(tauPrime * tauPrime + Math.cos(lambda - lambda0) ** 2));
  let xi = xiPrime;
  let eta = etaPrime;
  for (let j = 1; j <= 6; j++) {
    xi -= beta[j - 1] * Math.sin(2 * j * xiPrime) * Math.cosh(2 * j * etaPrime);
    eta -= beta[j - 1] * Math.cos(2 * j * xiPrime) * Math.sinh(2 * j * etaPrime);
  }
  return { x: K0 * Abar * (xi - meridionalArc(phi0) / Abar), y: K0 * Abar * eta };
}

function inversePlaneRectangular(x, y, zone) {
  const origin = zoneDefinitions[zone];
  if (!origin) throw new Error('座標系を選択してください。');
  const phi0 = origin.lat * RAD;
  const lambda0 = origin.lon * RAD;
  // Newton iteration against the exact Krüger forward projection, using a local numerical Jacobian.
  let phi = phi0 + x / (A * K0);
  let lambda = lambda0 + y / (A * K0 * Math.cos(phi0));
  for (let i = 0; i < 12; i++) {
    const p = forwardTM(phi, lambda, phi0, lambda0);
    const dx = x - p.x;
    const dy = y - p.y;
    if (Math.hypot(dx, dy) < 1e-7) break;
    const h = 1e-7;
    const pPhi = forwardTM(phi + h, lambda, phi0, lambda0);
    const pLam = forwardTM(phi, lambda + h, phi0, lambda0);
    const a = (pPhi.x - p.x) / h;
    const b = (pLam.x - p.x) / h;
    const c = (pPhi.y - p.y) / h;
    const d = (pLam.y - p.y) / h;
    const det = a * d - b * c;
    if (!Number.isFinite(det) || Math.abs(det) < 1e-8) throw new Error('入力座標を変換できません。X・Yと座標系を確認してください。');
    phi += (dx * d - b * dy) / det;
    lambda += (a * dy - dx * c) / det;
  }
  const lat = phi / RAD;
  const lon = lambda / RAD;
  if (!Number.isFinite(lat) || !Number.isFinite(lon) || Math.abs(lat) > 90 || Math.abs(lon) > 180) {
    throw new Error('変換結果が有効な範囲外です。入力値と座標系を確認してください。');
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
    if (dem !== null && dem >= -9999) setText('dem', `${fmt(dem)}${demResult.value.hsrc ? `（${demResult.value.hsrc}）` : ''}`);
    else setText('dem', '取得できませんでした');

    const failures = [geoid === null && 'ジオイド高', dem === null && '地理院標高'].filter(Boolean);
    setText('api-note', failures.length ? `${failures.join('・')}を取得できませんでした。通信環境を確認するか、国土地理院のサイトでご確認ください。` : 'ジオイド高は国土地理院API、地理院標高はDEMから取得しました。地形データの解像度や取得状況により値が異なる場合があります。');
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
    } catch { button.textContent = '選択してコピー'; }
  });
});
