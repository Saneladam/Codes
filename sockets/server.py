#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Mon 01. Jun 2026
#
# Purpose:      Listen for connections
# =============================================================================

# %% Import
import socket
import threading

# %% Constant
HOST = "0.0.0.0"
PORT = 9999

# %% Server
server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.bind((HOST, PORT))
server.listen()
clients = []  # (socket, name)


# %% Functions
def broadcast(msg):
    dead = []
    for client, _ in clients:
        try:
            client.send(msg)
        except:
            dead.append((client, _))
    for d in dead:
        clients.remove(d)


def handle_client(client, name):
    try:
        while True:
            data = client.recv(1024)
            if not data:
                break

            msg = f"{name}: {data.decode()}".encode()
            broadcast(msg)

    finally:
        print(f"{name} disconnected")
        clients.remove((client, name))
        client.close()


# %% Main
print("Listening...")
while True:
    client, addr = server.accept()
    print(f"Connected: {addr}")
    client.send(b"NAME?\n")
    name = client.recv(1024).decode().strip()
    clients.append((client, name))
    threading.Thread(target=handle_client, args=(client, name), daemon=True).start()
