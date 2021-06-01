from random import SystemRandom
from tqdm import trange, tqdm


class Dude:
    def __init__(self, curve, point, point_order):
        self.curve = curve
        self.point = point
        self.point_order = point_order


class FairDude(Dude):
    def __init__(self, curve, point, point_order):
        super().__init__(curve, point, point_order)
        self.cryptogen = SystemRandom()
        self.private_key = self.cryptogen.randint(2, self.point_order - 2)
        self.public_key = curve.mul(self.private_key, self.point)
        self.shared_key = None


class Alice(FairDude):
    def __init__(self, curve, point, point_order):
        super().__init__(curve, point, point_order)
        self.public_key_bob = None

    def send_public_key(self, bob, eva):
        bob.public_key_alice = self.public_key
        eva.public_key_alice = self.public_key

    def make_shared_key(self):
        self.shared_key = self.curve.mul(self.private_key, self.public_key_bob)

    def send_message(self, bob, eva):
        message = self.cryptogen.randint(2, self.curve.galois_field.p)
        encrypted_message = self.curve.galois_field.mul(message, self.shared_key.x)
        bob.encrypted_message_alice = encrypted_message
        eva.encrypted_message_alice = encrypted_message
        print('Alice sent message:', message, '(' + str(encrypted_message) + ' encrypted)')


class Bob(FairDude):
    def __init__(self, curve, point, point_order):
        super().__init__(curve, point, point_order)
        self.public_key_alice = None
        self.encrypted_message_alice = None

    def send_public_key(self, alice, eva):
        alice.public_key_bob = self.public_key
        eva.public_key_bob = self.public_key

    def make_shared_key(self):
        self.shared_key = self.curve.mul(self.private_key, self.public_key_alice)

    def decrypt_message(self):
        x_inverse = self.curve.galois_field.mul_inverse(self.shared_key.x)
        message = self.curve.galois_field.mul(self.encrypted_message_alice, x_inverse)
        print('Bob received message:', self.encrypted_message_alice, '(' + str(message) + ' decrypted)')


class Eva(Dude):
    def __init__(self, curve, point, point_order):
        super().__init__(curve, point, point_order)
        self.public_key_alice = None
        self.public_key_bob = None
        self.encrypted_message_alice = None
        self.shared_key = None
        self.private_key_alice = None

    def decrypt_message(self):
        potential_public_key_alice = self.point
        private_key_alice = 1

        for private_key_alice in trange(2, self.point_order, leave=False, desc='XAKNH7'):
            potential_public_key_alice = self.curve.add(potential_public_key_alice, self.point)
            if potential_public_key_alice == self.public_key_alice:
                break
        self.private_key_alice = private_key_alice

        self.shared_key = self.curve.mul(self.private_key_alice, self.public_key_bob)
        x_inverse = self.curve.galois_field.mul_inverse(self.shared_key.x)
        message = self.curve.galois_field.mul(self.encrypted_message_alice, x_inverse)
        print('Eva received message:', self.encrypted_message_alice, '(' + str(message) + ' decrypted)')
