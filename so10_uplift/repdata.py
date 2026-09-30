"""
Loader for the stored representation data, so later scripts need not
regenerate high levels.

    from so10_uplift.repdata import RepData
    D = RepData()
    D.max_level
    D.r8(n)            # SO(8) content of level n   (Sage D4 character)
    D.r9(n)            # SO(9) content r^9_n        (Sage B4 character)
    D.R9(n)            # cumulative sum_{i<=n} r^9_i (Sage B4 character)
    D.R10(n)           # canonical SO(10) uplift     (Sage D5 character)
    D.R10_extremal(n)  # {"min_height"|"max_l5"|"min_l5": D5 character}
    D.R10_all(n)       # list of all uplifts (D5 characters), n <= 10
    D.info(n)          # raw JSON record for level n of R10.json

Conjugation classes: uplifts are stored with l5 >= 0 representatives (Dynkin
a5 >= a4).  Replacing any lambda by its conjugate gives another uplift with
the same SO(9) restriction.  Raw JSON entries (ambient2 = doubled
orthogonal weight, dynkin, dim, mult) are available via the *_entries methods.
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


class RepData:
    def __init__(self, data_dir=DATA):
        with open(os.path.join(data_dir, "r9.json")) as f:
            self._r9 = json.load(f)
        with open(os.path.join(data_dir, "R10.json")) as f:
            self._R10 = json.load(f)
        self._all_path = os.path.join(data_dir, "R10_all_uplifts.json")
        self._all = None
        self.max_level = min(self._r9["meta"]["max_level"], self._R10["meta"]["max_level"])
        from sage.all import WeylCharacterRing
        self.D4 = WeylCharacterRing("D4", style="coroots")
        self.B4 = WeylCharacterRing("B4", style="coroots")
        self.D5 = WeylCharacterRing("D5", style="coroots")

    # raw entries
    def r8_entries(self, n):
        return self._r9["levels"][str(n)]["r8_D4"]

    def r9_entries(self, n):
        return self._r9["levels"][str(n)]["r9_B4"]

    def R9_entries(self, n):
        return self._R10["levels"][str(n)]["R9_B4"]

    def R10_entries(self, n):
        return self._R10["levels"][str(n)]["canonical_R10_D5"]

    def info(self, n):
        return self._R10["levels"][str(n)]

    # Sage characters
    @staticmethod
    def _char(ring, entries):
        c = ring(0)
        for e in entries:
            c += e["mult"] * ring(*e["dynkin"])
        return c

    def r8(self, n):
        return self._char(self.D4, self.r8_entries(n))

    def r9(self, n):
        return self._char(self.B4, self.r9_entries(n))

    def R9(self, n):
        return self._char(self.B4, self.R9_entries(n))

    def R10(self, n):
        return self._char(self.D5, self.R10_entries(n))

    def R10_extremal(self, n):
        return {k: self._char(self.D5, v)
                for k, v in self.info(n)["extremal_R10_D5"].items()}

    def R10_all(self, n):
        if self._all is None:
            with open(self._all_path) as f:
                self._all = json.load(f)["levels"]
        rec = self._all[str(n)]
        out = []
        for vec in rec["uplifts"]:
            c = self.D5(0)
            for k, dyn in zip(vec, rec["classes_dynkin"]):
                if k:
                    c += k * self.D5(*dyn)
            out.append(c)
        return out
