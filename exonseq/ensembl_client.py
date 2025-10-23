import time
from typing import Any, Dict, List, Optional

import requests


class EnsemblClient:
	def __init__(self, base_url: str = "https://rest.ensembl.org", timeout_s: int = 30, max_retries: int = 5):
		self.base_url = base_url.rstrip("/")
		self.timeout_s = timeout_s
		self.max_retries = max_retries
		self.session = requests.Session()
		self.session.headers.update({
			"Content-Type": "application/json",
			"Accept": "application/json",
			"User-Agent": "dystroSeq/1.0 (+https://example.org)"
		})

	def _request(self, method: str, path: str, params: Optional[Dict[str, Any]] = None, accept: str = "application/json") -> Any:
		url = f"{self.base_url}{path}"
		backoff = 1.0
		for attempt in range(1, self.max_retries + 1):
			headers = {**self.session.headers, "Accept": accept}
			resp = self.session.request(method, url, params=params, timeout=self.timeout_s, headers=headers)
			if resp.status_code == 429:
				retry_after = resp.headers.get("Retry-After")
				wait = float(retry_after) if retry_after else backoff
				time.sleep(wait)
				backoff = min(backoff * 2, 10.0)
				continue
			if 200 <= resp.status_code < 300:
				if accept == "application/json":
					return resp.json()
				return resp.text
			raise RuntimeError(f"Ensembl API error {resp.status_code}: {resp.text}")
		raise RuntimeError("Exhausted retries")

	def lookup_gene_by_symbol(self, species: str, symbol: str) -> Dict[str, Any]:
		return self._request("GET", f"/lookup/symbol/{species}/{symbol}", params={"expand": 1}, accept="application/json")

	def lookup_id(self, ensembl_id: str, expand: bool = True) -> Dict[str, Any]:
		params = {"expand": 1} if expand else None
		return self._request("GET", f"/lookup/id/{ensembl_id}", params=params, accept="application/json")

	def gene_transcripts(self, gene_id: str) -> List[Dict[str, Any]]:
		obj = self.lookup_id(gene_id, expand=True)
		transcripts = obj.get("Transcript") or obj.get("transcripts") or []
		return transcripts

	def sequence_region(self, species: str, region: str, coord_system_version: Optional[str] = None) -> str:
		params: Dict[str, Any] = {}
		if coord_system_version:
			params["coord_system_version"] = coord_system_version
		result = self._request("GET", f"/sequence/region/{species}/{region}", params=params, accept="text/plain")
		# Handle case where API returns JSON instead of plain text
		if result.startswith('{'):
			import json
			data = json.loads(result)
			return data.get("seq", "")
		return result

	def sequence_id(self, ensembl_id: str, seq_type: str) -> str:
		# seq_type in {"cdna","cds","protein"}
		result = self._request("GET", f"/sequence/id/{ensembl_id}", params={"type": seq_type}, accept="text/plain")
		# Handle case where API returns JSON instead of plain text
		if result.startswith('{'):
			import json
			data = json.loads(result)
			return data.get("seq", "")
		return result
