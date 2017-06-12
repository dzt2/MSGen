#include "bitseq.h"
#include "cerror.h"
#include <queue>
#include "ctext.h"

// BitSeq 
BitSeq::BitSeq(const BitSeq & seq) : bit_num(seq.bit_num), length(seq.length) {
	bytes = nullptr;
	if (seq.length > 0) {
		bytes = new byte[seq.length];
		for (int i = 0; i < seq.length; i++) {
			bytes[i] = seq.bytes[i];
		}
	}
}
BitSeq::BitSeq(BitSeq::size_t bitnum) : bit_num(bitnum) {
	length = bitnum / 8;
	if (bitnum % 8 != 0) length++;

	bytes = nullptr;
	if (length > 0) {
		bytes = new byte[length];
		for (int i = 0; i < length; i++) {
			bytes[i] = 0;
		}
	}
}
BitSeq::~BitSeq() {
	if (bytes != nullptr)
		delete bytes;
}
BitSeq::size_t BitSeq::bit_number() const { return bit_num; }
bit BitSeq::get_bit(BitSeq::size_t index) const {
	if (index >= bit_num)
		throw "Invalid index: ", index, " ( limits = ", bit_num, " )";
	else {
		byte bk = bytes[index / 8];
		return (bk & BIT_LOC[index % 8]) != 0;
	}
}
void BitSeq::set_bit(BitSeq::size_t index, bit val) {
	if (index >= bit_num)
		throw "Invalid index: ", index, " ( limits = ", bit_num, " )";
	else {
		byte bk = bytes[index / 8];
		bit oval = ((bk & BIT_LOC[index % 8]) != 0);
		if (oval != val)
			bytes[index / 8] = bytes[index / 8] ^ BIT_LOC[index % 8];
	}
}
std::string BitSeq::to_string() const {
	std::string str;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < 8 && i * 8 + j < bit_num; j++) {
			if ((bytes[i] & BIT_LOC[j]) != 0)
				str += '1';
			else str += '0';
		}
	}
	return str;
}
BitSeq BitSeq::subseq(BitSeq::size_t start, BitSeq::size_t end) const {
	if (start > bit_num)
		throw "Invalid index: ", start, " ( limits = ", bit_num, " )";
	else if (end > bit_num)
		"Invalid index: ", end, " ( limits = ", bit_num, " )";
	else if (start > end)
		"Invalid index: ", start, "; ", end;
	else {
		BitSeq seq(end - start); int k = 0;
		while (start < end) {
			bit obit = this->get_bit(start++);
			seq.set_bit(k++, obit);
		}
		return seq;
	}
}
bool BitSeq::subsume(const BitSeq & y) const {
	int n = (length < y.length) ? length : y.length;
	for (int i = 0; i < n; i++) {
		byte b1 = bytes[i];
		byte b2 = y.bytes[i];
		if (b1 != (b1 & b2)) return false;
	}
	return true;
}
bool BitSeq::equals(const BitSeq & seq) const {
	if (bit_num != seq.bit_num) return false;
	else {
		for (int i = 0; i < length; i++)
			if (bytes[i] != seq.bytes[i])
				return false;
		return true;
	}
}
int BitSeq::byte_number() const { return length; }
byte * BitSeq::get_bytes() const { return bytes; }
void BitSeq::set_bytes(const byte * new_bytes, size_t size) {
	size_t n = (length < size) ? length : size;
	for (int i = 0; i < n; i++)
		bytes[i] = new_bytes[i];
}
void BitSeq::clear_bytes() {
	for (int i = 0; i < length; i++)
		bytes[i] = 0;
}
void BitSeq::assign(const BitSeq & seq) {
	int bits = (bit_num < seq.bit_num) ? bit_num : seq.bit_num;
	for (int i = 0; i < bits; i++)
		this->set_bit(i, seq.get_bit(i));
}
void BitSeq::conjunct(const BitSeq & y) {
	if (this->bit_num != y.bit_num) {
		CError error(CErrorType::InvalidArguments, "BitSeq::conjunct", "Lenght is not matched: ("
			+ std::to_string(bit_num) + " <--> " + std::to_string(y.bit_num) + ")");
		CErrorConsumer::consume(error);
	}
	else {
		for (int i = 0; i < length; i++) 
			bytes[i] = (bytes[i] & (y.bytes)[i]);
	}
}
void BitSeq::disjunct(const BitSeq & y) {
	if (this->bit_num != y.bit_num) {
		CError error(CErrorType::InvalidArguments, "BitSeq::conjunct", "Lenght is not matched: ("
			+ std::to_string(bit_num) + " <--> " + std::to_string(y.bit_num) + ")");
		CErrorConsumer::consume(error);
	}
	else {
		for (int i = 0; i < length; i++)
			bytes[i] = (bytes[i] | (y.bytes)[i]);
	}
}
void BitSeq::bit_not() {
	for (int i = 0; i < length; i++)
		bytes[i] = ~(bytes[i]);
}
bool BitSeq::all_zeros() const {
	for (int i = 0; i < length; i++) {
		if (bytes[i] != '\0')
			return false;
	}
	return true;
}
BitSeq::size_t BitSeq::degree() const {
	size_t ones = 0, k = 0;
	while (k < bit_num) {
		if (get_bit(k++) == BIT_1)
			ones++;
	}
	return ones;
}

// BitTrie
BitTrie::BitTrie(BitSeq::size_t bias_index, const BitSeq & partial_key)
	: bias(bias_index), left(nullptr), right(nullptr), parent(nullptr), data(nullptr) {
	key = new BitSeq(partial_key);
}
BitTrie::~BitTrie() { delete key; }
BitSeq::size_t BitTrie::get_bias() const { return bias; }
const BitSeq & BitTrie::get_key() const { return *key; }
BitTrie * BitTrie::get_left() const { return left; }
BitTrie * BitTrie::get_right() const { return right; }
BitTrie * BitTrie::get_parent() const { return parent; }
void BitTrie::set_left(BitTrie * left_child) {
	if (left != nullptr)
		left->parent = nullptr;

	left = left_child;

	if (left != nullptr) {
		if (left->parent != nullptr) {
			if (left == left->parent->left) {
				left->parent->left = nullptr;
			}
			else {
				left->parent->right = nullptr;
			}
		}
		left->parent = this;
	}
}
void BitTrie::set_right(BitTrie * right_child) {
	if (right != nullptr)
		right->parent = nullptr;

	right = right_child;

	if (right != nullptr) {
		if (right->parent != nullptr) {
			if (right == right->parent->left)
				right->parent->left = nullptr;
			else
				right->parent->right = nullptr;
		}
		right->parent = this;
	}
}
bool BitTrie::is_leaf() const { return (left == nullptr) && (right == nullptr); }
void * BitTrie::get_data() const { return data; }
void BitTrie::set_data(void * value) { data = value; }

// BitTrieTree
BitTrieTree::BitTrieTree() : root(nullptr) {}
BitTrieTree::~BitTrieTree() {
	std::queue<BitTrie *> queue;
	if (root != nullptr) {
		queue.push(root);
		while (!queue.empty()) {
			BitTrie * node = queue.front();
			if (node->get_left() != nullptr) queue.push(node->get_left());
			if (node->get_right() != nullptr) queue.push(node->get_right());
			queue.pop(); delete node;
		}
	}
}
BitTrie * BitTrieTree::get_root() const { return root; }
BitTrie * BitTrieTree::maximum_prefix_match(const BitSeq & seq, BitSeq::size_t & index) const {
	// initialization 
	index = 0; BitTrie * node = root, *next;
	int i, seql, keyl; bit seqi, keyi;

	// match from 0 to specific node
	seql = seq.bit_number();
	while (node != nullptr) {
		/* match the bits in current node.key */
		keyl = (node->get_key()).bit_number();
		for (i = 0; i < keyl && index < seql; i++, index++) {
			seqi = seq.get_bit(index);
			keyi = (node->get_key()).get_bit(i);
			if (seqi != keyi) { break; }
		}

		/* not all-matched for this node or all-matched */
		if (i < keyl || index >= seql) break;
		/* all-matched for this node but not completed, to the next level */
		else {
			seqi = seq.get_bit(index++);
			if (seqi) next = node->get_right();
			else next = node->get_left();

			if (next == nullptr) break;
			else node = next;
		}
	} /* end while */

	  // return 
	return node;
}
BitTrie * BitTrieTree::get_leaf(const BitSeq & seq) const {
	BitSeq::size_t index = 0;
	BitTrie * leaf = this->maximum_prefix_match(seq, index);
	if (leaf != nullptr && index >= seq.bit_number())
		return leaf;
	else return nullptr;
}
BitTrie * BitTrieTree::insert_vector(const BitSeq & seq) {
	if (root == nullptr) {
		root = new BitTrie(0, seq);
		return root;
	}
	else {
		BitSeq::size_t index, bit_num = seq.bit_number();
		BitTrie * node = maximum_prefix_match(seq, index);
		if (index >= bit_num) return node;

		// get N, P, C1, C2
		BitTrie * parent = node->get_parent();
		BitTrie * lchild = node->get_left();
		BitTrie * rchild = node->get_right();

		// create prev, ori_post, new_post
		const BitSeq & key = node->get_key();
		BitSeq::size_t bias = index - (node->get_bias());
		BitSeq prev = key.subseq(0, bias);
		BitSeq ori_post = key.subseq(bias + 1, key.bit_number());
		BitSeq new_post = seq.subseq(index + 1, seq.bit_number());
		bit branch = seq.get_bit(index);

		// create N1, N2, L
		BitTrie * prev_node = new BitTrie(node->get_bias(), prev);
		BitTrie * post_node = new BitTrie(index + 1, ori_post);
		BitTrie * new_leaf = new BitTrie(index + 1, new_post);

		// connect N1 to N2 and L
		if (branch) {
			prev_node->set_left(post_node);
			prev_node->set_right(new_leaf);
		}
		else {
			prev_node->set_left(new_leaf);
			prev_node->set_right(post_node);
		}

		// connect N1 to P
		if (parent != nullptr) {
			if (parent->left == node)
				parent->set_left(prev_node);
			else
				parent->set_right(prev_node);
		}

		// connect N2 to C1 and C2
		post_node->set_left(lchild);
		post_node->set_right(rchild);

		// maintain the data
		if (node->data != nullptr) {
			post_node->data = node->data;
		}

		// remove node and update tree
		if (node == root)
			root = prev_node;
		delete node;

		// return 
		return new_leaf;
	}
}

void printTrieTree(const BitTrie & node, std::ostream & out, const std::string & prefix) {
	out << prefix << "\"" << node.get_key().to_string() << "\"\t(" << node.get_bias() << ")\n";
	if (!node.is_leaf()) {
		printTrieTree(*(node.get_left()), out, prefix + "|-- ");
		printTrieTree(*(node.get_right()), out, prefix + "|-- ");
	}
}