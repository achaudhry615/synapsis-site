"""JSONL run logs and numpy-free-readable snapshots."""
import json, os, time, hashlib, subprocess


def git_sha(default="unknown"):
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return default


class RunLog:
    def __init__(self, path, meta):
        self.path = path
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        self.f = open(path, "w")
        self.write({"record": "meta", "created": time.time(),
                    "git_sha": git_sha(), **meta})

    def write(self, rec):
        self.f.write(json.dumps(rec, default=float) + "\n")
        self.f.flush()

    def close(self):
        self.f.close()
        return sha256_file(self.path)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def read_jsonl(path):
    with open(path) as f:
        return [json.loads(l) for l in f if l.strip()]
