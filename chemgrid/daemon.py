import socketserver
import typer

class MyTCPHandler(socketserver.BaseRequestHandler):

    def handle(self):
        self.data = self.request.recv(1024).strip()
        print("{} wrote:".format(self.client_address[0]))
        print(self.data)

    @staticmethod
    def start(port: int):
        with socketserver.TCPServer(('localhost', port), MyTCPHandler) as server:
            server.serve_forever()

if __name__ == "__main__":
    MyTCPHandler.start(1555)
