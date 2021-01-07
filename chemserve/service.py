from pathlib import Path
from typing import Union, Optional, Callable, Mapping, Any, Sequence

from serviceit import ServiceClient, ServiceServer
from serviceit.server import Json

from chemserve.models import BaseChem, Concrete
from chemserve.payload import Payload, ConcretePayload


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
    def simple_server(
        cls, fn: Callable[[Concrete, Mapping[str, Any]], Any], port: int, return_hash: bool = False
    ) -> Server:
        def receiver(go: Json):
            payload = ConcretePayload.decode(go)
            values = {}
            for compound, item in zip(go["compounds"], payload.compounds):
                value = fn(item, payload.params)
                if return_hash:
                    base = BaseChem(compound["seq"], compound["name"], compound["key"])
                    values[hex(hash(base))] = value
                else:
                    values[compound["seq"]] = value
            return values

        return Server(receiver, port)

    @classmethod
    def client(cls, port: int) -> Client:
        return Client(port)

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


__all__ = ["Client", "Server", "ChemServe"]
