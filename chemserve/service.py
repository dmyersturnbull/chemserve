from __future__ import annotations
from pathlib import Path
from typing import Union, Optional, Callable, Mapping, Any, Sequence

from serviceit import ServiceClient, ServiceServer
from serviceit.server import Json

from chemserve.models import BaseChem, Concrete
from chemserve.payload import Payload, ConcretePayload


def say_compound(concrete: Concrete, params: Mapping[str, Any]) -> Payload:
    print(concrete.name + " " + concrete.inchikey)
    return Payload.empty()


def say_compounds(payload: ConcretePayload, params: Mapping[str, Any]) -> Payload:
    for concrete in payload.compounds:
        print(concrete.name + " " + concrete.inchikey)
    return Payload.empty()


class Client(ServiceClient):
    def send(self, payload: Union[Payload, Json]):
        super().send(payload)

    def send_one(self, chem: BaseChem, params=None):
        payload = Payload([chem], params).encode()
        super().send(payload)

    def send_multiple(self, chems: Sequence[BaseChem], params=None):
        payload = Payload(chems, params).encode()
        super().send(payload)


class Server(ServiceServer):
    def client(self) -> ServiceClient:
        return Client(self.port)


class DrawingServer(Server):
    @classmethod
    def draw(cls, payload: ConcretePayload) -> Json:
        from chemserve.graphics import ChemGraphicsKit

        chems = payload.compounds
        rows, cols = payload.bool_param("rows"), payload.bool_param("cols")
        path = Path(payload.str_param("path"))
        grid = ChemGraphicsKit().draw_grid(chems, rows, cols)
        grid.save(str(path))
        return Json({hex(hash(payload)): "success"})


class SimpleServerBuilder:
    def __init__(self, name: str = "server") -> None:
        self._name = name
        self._task_map = {}
        self._return_hash = False

    def map(
        self,
        name: str = "main",
        fn: Callable[[Concrete, Mapping[str, Any]], Any] = say_compound,
    ) -> SimpleServerBuilder:
        """
        Adds a task specified by some name.
        The name corresponds to a payload param called 'task' (aliased in ``Payload.task``).
        """
        self._task_map[name] = fn
        return self

    def return_by_hash(self) -> SimpleServerBuilder:
        """
        If True, the server returns a dict mapping from the hash of the compound
        to the return value of the function.
        Otherwise, the dict keys are the compound ``seq`` values.
        """
        self._return_hash = True
        return self

    def build(self, port: int) -> Server:
        if len(self._task_map) == 0:
            raise ValueError(f"No tasks added to {self._name}")

        def receiver(go: Json):
            payload = ConcretePayload.decode(go)
            values = {}
            fn = self._task_map[payload.task]
            for compound, item in zip(go["compounds"], payload.compounds):
                value = fn(item, payload.params)
                if self._return_hash:
                    base = BaseChem(compound["seq"], compound["name"], compound["key"])
                    values[hex(hash(base))] = value
                else:
                    values[compound["seq"]] = value
            return values

        return Server(receiver, port)


class FullServerBuilder:
    def __init__(self, name: str = "server") -> None:
        self._name = name
        self._task_map = {}

    def map(
        self, name: str = "main", fn: Callable[[ConcretePayload], Any] = say_compounds
    ) -> FullServerBuilder:
        """
        Adds a task specified by some name.
        The name corresponds to a payload param called 'task' (aliased in ``Payload.task``).
        """
        self._task_map[name] = fn
        return self

    def build(self, port: int) -> Server:
        if len(self._task_map) == 0:
            raise ValueError(f"No tasks added to {self._name}")

        def receiver(go: Json):
            payload = ConcretePayload.decode(go)
            fn = self._task_map[payload.task]
            return fn(payload)

        return Server(receiver, port)


class ChemServe:
    @classmethod
    def drawing_server(cls, port: int) -> DrawingServer:
        def receiver(go: Json):
            payload = ConcretePayload.decode(go)
            return DrawingServer.draw(payload)

        return DrawingServer(receiver, port=port)

    @classmethod
    def server(cls, fn: Callable[[ConcretePayload], Optional[Json]], port: int) -> Server:
        def receiver(go: Json):
            payload = ConcretePayload.decode(go)
            return fn(payload)

        return Server(receiver, port)

    @classmethod
    def per_compound(cls, name: str) -> SimpleServerBuilder:
        return SimpleServerBuilder(name)

    @classmethod
    def many_compounds(cls, name: str) -> FullServerBuilder:
        return FullServerBuilder(name)

    @classmethod
    def client(cls, port: int) -> Client:
        return Client(port)

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


__all__ = ["Client", "Server", "ChemServe"]
