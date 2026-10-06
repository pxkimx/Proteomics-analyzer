import json
import os
import signal
import subprocess
import sys
import time
import urllib.request

import pytest

PORT = 8776
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = f"http://127.0.0.1:{PORT}"


def call(path, body=None, raw=None):
    data = raw if raw is not None else (None if body is None else json.dumps(body).encode())
    req = urllib.request.Request(BASE + path, data=data, method="POST" if data is not None else "GET")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def wait_up():
    for _ in range(60):
        try:
            return call("/api/state")[0] == 200
        except Exception:
            time.sleep(0.1)
    return False


@pytest.fixture()
def server():
    env = dict(os.environ, RD_PARENT_PID=str(os.getpid()))
    p = subprocess.Popen([sys.executable, "-m", "recode_detector.server", "--port", str(PORT), "--no-browser"],
                         cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert wait_up()
    yield p
    if p.poll() is None:
        p.send_signal(signal.SIGTERM)
        p.wait(5)


def test_full_workflow_over_http(server):
    code, body = call("/api/example", {})
    st = json.loads(body)
    assert code == 200 and st["loaded"] and len(st["samples"]) == 8 and st["kinds"] == ["LFQ intensity"]
    code, body = call("/api/compare", {"a": "Treated", "b": "Control"})
    assert code == 400 and "processing" in json.loads(body)["error"]
    code, body = call("/api/analyze", {"params": {"normalization": "median"}})
    r = json.loads(body)
    assert code == 200 and r["groups"] == ["Control", "Treated"] and r["qc"]["n_proteins"] > 2000
    code, body = call("/api/compare", {"a": "Treated", "b": "Control", "fdr": 0.05, "lfc": 1.0})
    c = json.loads(body)
    assert code == 200 and c["counts"]["up"] + c["counts"]["down"] > 100 and len(c["heatmap"]["z"]) > 0
    hit = next(r for r in c["rows"] if r["call"] != "ns")
    code, body = call("/api/protein?id=" + hit["protein"])
    assert code == 200 and len(json.loads(body)["values"]) == 8
    code, body = call("/api/download/results")
    assert code == 200 and body.startswith(b"protein,log2fc")
    code, body = call("/api/download/settings")
    assert json.loads(body)["parameters"]["normalization"] == "median"
    gmt = "S1\tx\t" + "\t".join(r["gene"] for r in c["rows"][:12]) + "\n"
    code, body = call("/api/enrich", {"gmt": gmt, "direction": "both"})
    assert code == 200 and json.loads(body)["rows"][0]["set"] == "S1"
    # the page itself and path traversal
    assert call("/")[0] == 200 and b"Recode Detector" in call("/")[1]
    assert call("/../recode_detector/server.py")[0] == 404


def test_bad_uploads_give_clear_errors(server):
    code, body = call("/api/load?name=x.txt", raw=b"hello")
    assert code == 400 and "error" in json.loads(body)
    assert call("/api/state")[0] == 200            # still alive


def test_quits_when_window_closes(server):
    import http.client
    conn = http.client.HTTPConnection("127.0.0.1", PORT)
    conn.request("GET", "/api/window")
    resp = conn.getresponse(); resp.readline()
    time.sleep(8)
    assert server.poll() is None                   # window open: keeps running
    resp.close(); conn.close()                     # tab closed
    t0 = time.time()
    while server.poll() is None and time.time() - t0 < 15:
        time.sleep(0.3)
    assert server.poll() is not None and time.time() - t0 < 12
    with pytest.raises(Exception):
        call("/api/state")                         # nothing listening any more


def test_quits_when_launcher_dies():
    launcher = subprocess.Popen(["sleep", "60"])
    env = dict(os.environ, RD_PARENT_PID=str(launcher.pid))
    p = subprocess.Popen([sys.executable, "-m", "recode_detector.server", "--port", str(PORT), "--no-browser"],
                         cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        assert wait_up()
        launcher.kill(); launcher.wait()
        p.wait(timeout=5)
    finally:
        if p.poll() is None:
            p.kill()
