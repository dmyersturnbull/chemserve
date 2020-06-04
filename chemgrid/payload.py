from __future__ import annotations
from typing import Union, Sequence, Mapping, Optional

from serviceit.server import Json, Payload as _Payload

from chemgrid.base_chem import BaseChem
from chemgrid.models import Concrete, BaseChem


JsonDataType = Union[int, str, None, bool, Sequence[int], Sequence[str], Sequence[bool]]


class Payload:
    def __init__(self, items: Sequence[BaseChem], params: Optional[Mapping[str, JsonDataType]]):
        self._items = items
        self._params = {} if params is None else params

    @classmethod
    def new(cls, items: Sequence[BaseChem], params: Mapping[str, JsonDataType]) -> Payload:
        return Payload(items, params)

    def encode(self) -> Json:
        return _Payload({
            'compounds': [self._encode(i) for i in self._items],
            'params': self._params
        })

    @property
    def items(self) -> Sequence[BaseChem]:
        return self._items

    @property
    def params(self) -> Mapping[str, JsonDataType]:
        return self._params

    def str_param(self, key: str) -> str:
        item = self._params[key]
        if not isinstance(item, str):
            raise ValueError("Param {} (value {}) has type {}".format(key, item, type(item)))
        return item

    def int_param(self, key: str) -> int:
        item = self._params[key]
        if not isinstance(item, int):
            raise ValueError("Param {} (value {}) has type {}".format(key, item, type(item)))
        return item

    def bool_param(self, key: str) -> bool:
        item = self._params[key]
        if not isinstance(item, bool):
            raise ValueError("Param {} (value {}) has type {}".format(key, item, type(item)))
        return item

    @classmethod
    def decode(cls, payload: Json) -> ConcretePayload:
        items = [cls._decode(c) for c in payload['compounds']]
        params = payload.get('params', {})
        return ConcretePayload(items, params)

    def _encode(self, item: BaseChem):
        return {
            'seq': item.inchi_or_smiles,
            'name': item.name,
            'key': item.key
        }

    @classmethod
    def _decode(cls, dct: Mapping[str, str]) -> Concrete:
        return Concrete.of(dct['seq'], dct.get('name'), dct.get('key'))

    def __getitem__(self, item: Union[int, str]):
        if isinstance(item, int):
            return self._items[item]
        elif isinstance(item, str) and item=='params':
            return self._params[item]
        else:
            raise TypeError("{} not found".format(item))

    def __eq__(self, other):
        # this is slightly harsh
        if not isinstance(other, self.__class__):
            raise TypeError("Universal equality unsupported")
        return self.params == other.params and self.items == other.items

    def __hash__(self):
        chem_hashes = [hash(item) for item in self.items]
        param_hashes = [str(k)+'='+str(v) for k, v in self.params.items()]
        return hash(tuple(*chem_hashes, *param_hashes))

    def __repr__(self):
        return "{}(n={} : ({}) @ {})".format(self.__class__.__name__, len(self.items), ','.join(k+'='+str(v) for k, v in self.params), hex(id(self)))

    def __str__(self):
        return "{}(n={} : ({}))".format(self.__class__.__name__, len(self.items), ','.join(k+'='+str(v) for k, v in self.params))


class ConcretePayload(Payload):

    @property
    def items(self) -> Sequence[Concrete]:
        # noinspection PyTypeChecker
        return self._items



__all__ = ['Payload', 'ConcretePayload']
