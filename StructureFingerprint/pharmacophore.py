from typing import Iterator, Tuple, TYPE_CHECKING, Union
from .feature_queries import *

if TYPE_CHECKING:
    from CGRtools.containers import MoleculeContainer, CGRContainer, QueryContainer

queries = {
    'acc_allow': acceptor_queries,
    'acc_ban': acceptor_banned,
    'acid_allow': acidic_queries,
    'ar_allow': aromatic_queries,
    'base_allow': basic_queries,
    'base_ban': basic_banned,
    'donor_allow': donor_queries,
    'hal_allow': halogen_queries,
}


class Features:
    __slots__ = ()

    @staticmethod
    def matched_ids(mol, query):
        return {dct[1] for dct in query.get_mapping(mol)}

    def features(self, mol) -> Iterator[Tuple[int, int]]:
        id_sets = {}
        for name, qrs in queries.items():
            if len(qrs) == 1:
                id_sets[name] = self.matched_ids(mol, qrs[0])
                continue
            first, others = qrs[0], qrs[1:]
            first = self.matched_ids(mol, first)
            first.union(*(self.matched_ids(mol, q) for q in others))
            id_sets[name] = first

        bins = [
            id_sets['acc_allow'] - id_sets['acc_ban'],
            id_sets['acid_allow'],
            id_sets['ar_allow'],
            id_sets['base_allow'] - id_sets['base_ban'],
            id_sets['donor_allow'],
            id_sets['hal_allow'],
        ]

        for idx, _ in mol.atoms():
            yield idx, int(''.join('1' if idx in set_ else '0' for set_ in bins), base=2)


__all__ = ['Features']
