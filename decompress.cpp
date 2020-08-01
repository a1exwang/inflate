#include <cassert>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <functional>
#include <tuple>
#include <vector>
#include <list>
#include <unordered_map>

#include <coroutine>

class InvalidFormat :public std::runtime_error {
 public:
  explicit InvalidFormat(const std::string &s) :runtime_error(s) { }
};

int get_byte(std::istream &is) {
  auto c = is.get();
  if (c == EOF) {
    throw std::out_of_range("stream EOF");
  }
  return c;
}
void assert_byte(std::istream &is, uint8_t expected) {
  auto c = get_byte(is);
  if (c != expected) {
    throw InvalidFormat("byte should be " + std::to_string((int)expected) + " but '" + std::to_string((int)c));
  }
}

void skip(std::istream &is, size_t n) {
  is.seekg(n, std::ios::cur);
  if (!is) {
    throw std::out_of_range("stream EOF");
  }
}

uint32_t read_uint32(std::istream &is) {
  uint32_t result = 0;
  result |= (uint32_t)is.get();
  result |= (uint32_t)is.get() << 8u;
  result |= (uint32_t)is.get() << 16u;
  result |= (uint32_t)is.get() << 24u;
  if (!is) {
    throw std::out_of_range("stream EOF");
  }
  return result;
}

void write_uint32(std::ostream &os, uint32_t value) {
  os.put(value & 0xffu); value >>= 8u;
  os.put(value & 0xffu); value >>= 8u;
  os.put(value & 0xffu); value >>= 8u;
  os.put(value & 0xffu);
}

constexpr size_t length_table[][2] = {
    {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9}, {0, 10}, {1, 11}, {1, 13}, {1, 15}, {1, 17}, {2, 19}, {2, 23}, {2, 27}, {2, 31}, {3, 35}, {3, 43}, {3, 51}, {3, 59}, {4, 67}, {4, 83}, {4, 99}, {4, 115}, {5, 131}, {5, 163}, {5, 195}, {5, 227}, {0, 258}, {std::numeric_limits<size_t>::max(), 259/* sentry object */}
};

// code -> (extra_bits, start_offset)
constexpr size_t distance_table[][2] = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 5}, {1, 7}, {2, 9}, {2, 13}, {3, 17}, {3, 25}, {4, 33}, {4, 49}, {5, 65}, {5, 97}, {6, 129}, {6, 193}, {7, 257}, {7, 385}, {8, 513}, {8, 769}, {9, 1025}, {9, 1537}, {10, 2049}, {10, 3073}, {11, 4097}, {11, 6145}, {12, 8193}, {12, 12289}, {13, 16385}, {13, 24577},{std::numeric_limits<size_t>::max(),32769}, {0,0}
};

class BitStreamReader {
 public:
  explicit BitStreamReader(std::istream &is) :is_(is) { }
  void ensure_cached_byte() {
    if (cached_byte_size == 0) {
      cached_byte_size = 8;
      auto c = is_.get();
      offset_++;
      if (c == EOF) {
        throw InvalidFormat("stream EOF");
      }
      cached_byte = c;
    }
  }
  void discard_remaining_bits() {
    cached_byte_size = 0;
  }

  // n < 32
  // 11000101 -> 10100011
  uint32_t read_bit(size_t n) {
    assert(n <= 32);
    uint32_t result = 0;
    for (size_t i = 0; i < n; i++) {
      ensure_cached_byte();
      result <<= 1u;
      result |= ((uint32_t)cached_byte >> (8-cached_byte_size)) & 1u;
      cached_byte_size--;
    }
    return result;
  }

  uint32_t read_bit_le(size_t n) {
    assert(n <= 32);
    uint32_t result = 0;
    for (size_t i = 0; i < n; i++) {
      ensure_cached_byte();
      auto target_bit = ((uint32_t)cached_byte >> (8-cached_byte_size)) & 1u;
      result |= (target_bit << i);
      cached_byte_size--;
    }
    return result;
  }

  void copy(std::ostream &os, size_t size) {
    for (size_t i = 0; i < size; i++) {
      os.put(is_.get());
      offset_++;
      if (!is_) {
        throw InvalidFormat("Unexpected EOF when copying non-compressed data");
      }
    }
  }

  template <typename T>
  T read_le() {
    T value;
    for (int i = 0; i < sizeof(T); i++) {
      reinterpret_cast<uint8_t*>(&value)[i] = is_.get();
      offset_++;
      assert(is_);
    }
    return value;
  }

  std::string offset_string() const {
    std::stringstream ss;
    if (cached_byte_size == 0) {
      ss << (offset_ + 1) << "+0b";
    } else {
      ss << offset_ << "+" << (8 - cached_byte_size) << "b";
    }
    return ss.str();
  }

 private:
  uint8_t cached_byte = 0;
  size_t cached_byte_size = 0;
  std::istream& is_;
  size_t offset_ = std::numeric_limits<size_t>::max();
};

class BitStreamWriter {
 public:
  explicit BitStreamWriter(std::ostream &os) :os(os) { }

  void flush_all() {
    if (cached_byte_size != 0) {
      os.put(byte);
      cached_byte_size = 0;
      byte = 0;
    }
  }

  void write_le(uint32_t bits, size_t n) {
    assert(n <= 32);
    uint32_t result = 0;
    for (size_t i = 0; i < n; i++) {
      auto target_bit = bits & 1u;
      bits >>= 1u;
      byte |= target_bit << cached_byte_size;
      cached_byte_size++;
      if (cached_byte_size == 8) {
        os.put(byte);
        cached_byte_size = 0;
        byte = 0;
      }
    }
  }

  void write_be(uint32_t bits, size_t n) {
    assert(n <= 32);
    uint32_t result = 0;
    for (size_t i = 0; i < n; i++) {
      auto target_bit = (bits >> (n-i-1)) & 1u;
      byte |= target_bit << cached_byte_size;
      cached_byte_size++;
      if (cached_byte_size == 8) {
        os.put(byte);
        cached_byte_size = 0;
        byte = 0;
      }
    }
  }
 private:
  size_t cached_byte_size = 0;
  uint8_t byte = 0;
  std::ostream &os;
};

struct HuffmanTree {
  HuffmanTree() = default;
  explicit HuffmanTree(const std::vector<size_t> &code_lengths)
      :counts(65, 0), lengths(code_lengths), bits(code_lengths.size()) {

    std::vector<std::vector<int>> n_codes(code_lengths.size(), std::vector<int>(0));
    for (int code = 0; code < code_lengths.size(); code++) {
      counts[code_lengths[code]]++;
      n_codes[code_lengths[code]].push_back(code);
    }
    while (!counts.empty() && counts.back() == 0) {
      if (counts.back() == 0) {
        counts.pop_back();
      }
    }

    offsets.resize(counts.size() + 1);
    skip.resize(offsets.size());

    if (!counts.empty()) {
      counts[0] = 0;
      for (int i = 1; i <= counts.size(); i++) {
        offsets[i] = counts[i-1] + offsets[i-1];
      }

      uint64_t last_skip = 0;
      for (int i = 1; i < counts.size(); i++) {
        last_skip = (last_skip + counts[i-1]) << 1u;
        skip[i] = last_skip;
      }

      size_t prev = 0;
      for (int i = 1; i < n_codes.size(); i++) {
        for (auto code : n_codes[i]) {
          bits[code] = prev;
          prev++;
          symbols.push_back(code);
        }
        prev <<= 1;
      }

    }
  }

  size_t size() const {
    return lengths.size();
  }

  std::tuple<size_t, size_t> at(size_t value) const {
    return std::make_tuple(bits[value], lengths[value]);
  }

  // Read a value from stream decoded by 'tree'
  uint32_t read_value(BitStreamReader &reader) const {
    uint64_t value = 0;
    size_t bit_length = 0;
    while (true) {
      value |= reader.read_bit_le(1);
      bit_length++;
      if (bit_length > offsets.size()) {
        throw InvalidFormat("bit_length overflows to tree.offset.size() = '" + std::to_string(bit_length) + "'");
      }

      auto off = value - skip[bit_length];
      if (off < counts[bit_length]) {
        auto index = offsets[bit_length] + off;
        return symbols[index];
      }
      value <<= 1u;
    }
  }
//  void print(std::ostream &os) {
//    for (int i = 0; i < symbols.size(); i++) {
//      os << "symbols[" << i << "] = " << symbols[i] << std::endl;
//    }
//    for (auto offset : offsets) {
//      os << offset << " ";
//    }
//    std::cout << std::endl;
//    for (auto count : counts) {
//      os << count << " ";
//    }
//    os << std::endl;
//  }

  void write_value(BitStreamWriter &writer, size_t code) const {
    assert(lengths[code] > 0);
    writer.write_be(bits[code], lengths[code]);
  }

  std::vector<size_t> offsets;
  std::vector<size_t> counts;
  std::vector<size_t> symbols;
  std::vector<uint64_t> skip;

  std::vector<size_t> bits, lengths;
};


template <typename T>
class RingBuffer {
 public:
  explicit RingBuffer(size_t max_size) :data_(max_size) {
    assert(max_size > 1);
  }

  T pop_front() {
    // TODO
    assert(0);
  }

  T at(size_t offset) const {
    if (start_ + offset < data_.size()) {
      return data_[start_ + offset];
    } else {
      auto offset_from0 = offset - (data_.size() - start_);
      return data_[offset_from0];
    }
  }

  // Iterate each item in the ring buffer, cb(value, offset)
  void iterate(std::function<bool (T, size_t)> cb) {
    if (start_ <= end_) {
      for (size_t i = start_; i < end_; i++) {
        cb(data_[i], i - start_);
      }
    } else {
      for (size_t i = start_; i < data_.size(); i++) {
        cb(data_[i], i - start_);
      }
      for (size_t i = 0; i < end_; i++) {
        cb(data_[i], i + (data_.size() - start_));
      }
    }
  }

  size_t size() const {
    if (start_ <= end_) {
      return end_ - start_;
    } else {
      return end_ + data_.size() - start_;
    }
  }

  // push one element into end_, if the buffer is full, delete the first element.
  std::optional<T> push_back(T data) {
    if (next(end_) == start_) {
      // buffer is full
      auto ret = data_[start_];
      start_ = next(start_);
      end_ = next(end_);
      data_[end_] = data;
      return ret;
    } else {
      data_[end_] = data;
      end_ = next(end_);
      return {};
    }
  }

  // duplicate from previous distance'th position, with the size of "size", into end_
  std::vector<std::span<const T>> dup(size_t distance, size_t size) {
    auto start_offset = prev(end_, distance);
    auto j = start_offset;
    for (size_t i = 0; i < size; i++, j = next(j)) {
      push_back(data_[j]);
    }
    if (start_offset + size <= data_.size()) {
      return {{&data_[start_offset], size}};
    } else {
      auto remaining_size = size - (data_.size() - start_offset);
      return {{data_.begin()+start_offset, data_.end()}, {data_.begin(), data_.begin()+remaining_size}};
    }
  }

  // return the next index of offset
  size_t next(size_t offset, size_t distance = 1) {
    return (offset + distance) % data_.size();
  }

  // return the previous distance'th index of offset
  size_t prev(size_t offset, size_t distance = 1) {
    assert(distance <= data_.size());
    return (offset + data_.size() - distance) % data_.size();
  }

 private:
  std::vector<T> data_;

  // start_: the first valid data index.
  // end_: the first empty slot to insert new data.
  // When start == end, the buffer is empty rather than full.
  // The maximum size of the buffer is data_.size()-1.
  size_t start_ = 0;
  size_t end_ = 0;
};

class Inflate {
 public:
  explicit Inflate(std::istream& is) :reader_(is), ring_buffer_(65*1024) { }

  size_t decompress(std::ostream &os) {
    size_t total_size = 0;
    while (true) {
//      std::cerr << "start to decompress output_off=" << total_size
//                << ", input_off=" << reader_.offset_string() << std::endl;
      auto t0 = std::chrono::high_resolution_clock::now();
      auto [size, eof] = decompress_block(os);
      auto t1 = std::chrono::high_resolution_clock::now();
//      std::cerr << "Block decode time " << std::chrono::duration<double, std::milli>(t1 - t0).count() << "ms" << std::endl;
      total_size += size;
      if (eof) {
        break;
      }
    }
    return total_size;
  }

  std::tuple<size_t, bool> decompress_block(std::ostream &os) {
    auto bfinal = reader_.read_bit_le(1);
    auto btype = reader_.read_bit_le(2);

    size_t block_size = 0;
    if (btype == 0b01) {
      // fixed Huffman tree
      block_size = decompress_block_data(os, true, HuffmanTree(), HuffmanTree());
    } else if (btype == 0b10) {
      // dynamic Huffman tree
      auto [lit_tree, dist_tree] = read_dynamic_huffman_tree_header();
      block_size = decompress_block_data(os, false, lit_tree, dist_tree);
    } else if (btype == 0) {
      // non-compressed data
      reader_.discard_remaining_bits();

      auto len = reader_.read_le<uint16_t>();
      auto nlen = reader_.read_le<uint16_t>();
      if ((len ^ nlen) != (uint16_t)0xffff) {
        throw InvalidFormat("Data corruption: len(" + std::to_string(len) + ") != ~nlen(" + std::to_string(nlen) + ")");
      }

      reader_.copy(os, len);

      block_size = len;
    } else {
      throw InvalidFormat("Unexpected block type 0b11");
    }

    return {block_size, bfinal};
  }

 private:

  // Read dynamic huffman code length table decoded from 'huffman_tree'
  std::vector<size_t> generate_huffman_table(HuffmanTree &huffman_tree, const size_t lit_count) {
    std::vector<size_t> huffman_table;
    huffman_table.reserve(lit_count);
    int prev = -1;

    for (size_t lit = 0; lit < lit_count; ) {
      auto value = huffman_tree.read_value(reader_);
      /**
       0 - 15: Represent code lengths of 0 - 15
       16: Copy the previous code length 3 - 6 times.
           The next 2 bits indicate repeat length
                 (0 = 3, ... , 3 = 6)
              Example:  Codes 8, 16 (+2 bits 11),
                        16 (+2 bits 10) will expand to
                        12 code lengths of 8 (1 + 6 + 5)
       17: Repeat a code length of 0 for 3 - 10 times.
           (3 bits of length)
       18: Repeat a code length of 0 for 11 - 138 times
           (7 bits of length)
       */
      if (value == 16) {
        if (prev == -1) {
          throw InvalidFormat("Should not use 'copy' code at the beginning");
        }
        auto copy_previous_count = reader_.read_bit_le(2) + 3;
        for (int i = 0; i < copy_previous_count; i++) {
          huffman_table.push_back(prev);
        }
        lit += copy_previous_count;
      } else if (17 <= value && value <= 18) {
        size_t zero_count = 0;
        if (value == 17) {
          zero_count = reader_.read_bit_le(3) + 3;
        } else {
          zero_count = reader_.read_bit_le(7) + 11;
        }
        for (int i = 0; i < zero_count; i++) {
          huffman_table.push_back(0);
        }
        prev = 0;
        lit += zero_count;
      } else if (value < 16) {
        huffman_table.push_back(value);
        prev = value;
        lit++;
      } else {
        throw std::runtime_error("decode invalid value, should be [0, 18], but '" + std::to_string(value) + "'");
      }
    }
//    for (int i = 0; i < huffman_table.size(); i++) {
//      if (huffman_table[i] != 0) {
//        std::cout << std::hex << std::setw(2) << std::setfill('0') << i << "(" << (char)i << ") --> " << huffman_table[i] << std::endl;
//      }
//    }

    return huffman_table;
  }

  // Read dynamic Huffman tree from block header
  std::tuple<HuffmanTree, HuffmanTree> read_dynamic_huffman_tree_header() {
    const auto hlit = reader_.read_bit_le(5) + 257;
    const auto hdist = reader_.read_bit_le(5) + 1;
    const auto hclen = reader_.read_bit_le(4) + 4;

    std::vector<size_t> code_lengths(19, 0);

    int index[] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
    for (int i = 0; i < hclen; i++) {
      auto code_length = reader_.read_bit_le(3);
      code_lengths[index[i]] = code_length;
    }

    HuffmanTree code_length_huffman_tree(code_lengths);

    auto lit_table = generate_huffman_table(code_length_huffman_tree, hlit);
    auto dist_table = generate_huffman_table(code_length_huffman_tree, hdist);

    return std::make_tuple(HuffmanTree(lit_table), HuffmanTree(dist_table));
  }

 private:
  size_t decompress_block_data(
      std::ostream &os,
      bool fixed,
      const HuffmanTree &lit_tree,
      const HuffmanTree &dist_tree) {

    size_t block_offset = 0;
    while (true) {

      uint32_t code;
      if (fixed) {
        code = read_fixed_lit();
      } else {
        code = lit_tree.read_value(reader_);
      }

      if (code >= 286) {
        throw InvalidFormat("Invalid code length code " + std::to_string(code));
      } else if (code == 256) {
        // end of block
        break;
      } else if (code < 256) {
        // literal code
        os.put(code);
        ring_buffer_.push_back(code);
        block_offset++;
      } else {
        // length/distance code

        // first read length
        size_t length = 0, distance = 0;
        auto length_code = code - 257;
        {
          auto [extra_bits, start_off] = length_table[length_code];
          length = reader_.read_bit_le(extra_bits) + start_off;
        }

        // then distance
        uint32_t dist_code = 0;
        if (fixed) {
          dist_code = reader_.read_bit(5);
        } else {
          dist_code = dist_tree.read_value(reader_);
        }

        {
          auto [extra_bits, start_off] = distance_table[dist_code];
          distance = reader_.read_bit_le(extra_bits) + start_off;
        }

        if (distance == 0) {
          throw InvalidFormat("distance should not be 0");
        }

        auto buffers = ring_buffer_.dup(distance, length);
        for (auto buffer : buffers) {
          os.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
        }
        block_offset += length;
      }
    }

    auto block_size = block_offset;

    return block_size;
  }


  uint32_t read_fixed_lit() {
    auto prefix4 = reader_.read_bit(4);
    if (prefix4 < 0b0011) {
      // 256 - 279
      auto part3 = reader_.read_bit(3);
      if (prefix4 & 1u) {
        part3 += 8;
      }
      return part3 + 256;
    } else if (prefix4 >= 0b0011 && prefix4 < 0b1100) {
      auto prefix5 = (prefix4 << 1u) | reader_.read_bit(1);
      if (prefix5 >= 0b00110 && prefix5 < 0b11000) {
        // 0 - 143
        return ((prefix5 << 3u) | reader_.read_bit(3)) - 0x30;
      } else if (prefix5 >= 0b11000 && prefix5 < 0b11001) {
        // 280 - 287
        return reader_.read_bit(3) + 280;
      } else {
        // 144 - 255, code starts from 0b110010000
        return ((prefix5 << 4u) | reader_.read_bit(4)) - 0b110010000 + 144;
      }
    } else {
      assert(0);
    }
  }

 private:
  BitStreamReader reader_;
  RingBuffer<uint8_t> ring_buffer_;
};

// RFC1952
void decompress_stream(std::istream &is, std::ostream &os) {

  // ID1
  assert_byte(is, 0x1f);
  // ID2
  assert_byte(is, 0x8b);
  // CM: compression method
  assert_byte(is, 0x08);

  /**
   bit 0   FTEXT
   bit 1   FHCRC
   bit 2   FEXTRA
   bit 3   FNAME
   bit 4   FCOMMENT
   bit 5   reserved
   bit 6   reserved
   bit 7   reserved
   */
  auto flags = get_byte(is);
  if (flags != 0) {
    throw InvalidFormat("Non-zero flag not implemented");
  }

  // mtime
  uint32_t mtime = read_uint32(is);
  {
    std::time_t temp = mtime;
    std::tm* t = std::gmtime(&temp);
    std::stringstream ss; // or if you're going to print, just input directly into the output stream
    ss << std::put_time(t, "%Y-%m-%d %I:%M:%S %p");
    std::string mtime_string = ss.str();
  }

  // extra flag
  skip(is, 1);

  // OS
  skip(is, 1);

  Inflate inflate(is);
  size_t decompressed_size = inflate.decompress(os);

  auto crc32 = read_uint32(is);
  (void)crc32;
  // TODO: verify CRC32
//  std::cout << "crc32: " << crc32 << std::endl;
  auto original_size = read_uint32(is);
  if (decompressed_size != original_size) {
    throw InvalidFormat("data corruption, original_size != decompressed_size");
  }
}

struct Node {
  Node(Node *left, Node *right) :left(left), right(right), freq(left->freq + right->freq) { }
  Node(size_t code, size_t freq) :code(code), freq(freq) { }
  bool operator<(const Node &rhs) const {
    if (freq == rhs.freq) {
      return code < rhs.code;
    } else {
      return freq > rhs.freq;
    }
  }
  bool is_leaf() const {
    return left == nullptr && right == nullptr;
  }

  void dfs(std::vector<size_t> &result, size_t depth = 0) const {
    if (is_leaf()) {
      result[code] = depth;
    } else {
      left->dfs(result, depth+1);
      right->dfs(result, depth+1);
    }
  }

  size_t freq = 0;
  size_t code = std::numeric_limits<size_t>::max();
  Node *left = nullptr;
  Node *right = nullptr;
};

// <bits, nbits>
HuffmanTree create_huffman_tree(const std::vector<size_t> &freqs) {
  std::vector<Node> nodes_owner;
  std::vector<Node*> nodes;
  // reserve memory to prevent invalidating pointer
  nodes_owner.reserve(2*freqs.size());
  nodes.reserve(2*freqs.size());
  for (size_t i = 0; i < freqs.size(); i++) {
    if (freqs[i] > 0) {
      nodes_owner.emplace_back(i, freqs[i]);
      nodes.emplace_back(&nodes_owner.back());
    }
  }

  std::vector<size_t> lengths(freqs.size(), 0);
  if (nodes.empty()) {
    // TODO: not implemented
    assert(0);
  } else if (nodes.size() == 1) {
    // DEFLATE format does not support zero length tree
    lengths[nodes.front()->code] = 1;
  } else {
    auto compare = [](Node* lhs, Node *rhs) { return *lhs < *rhs; };
    std::make_heap(nodes.begin(), nodes.end(), compare);
    while (nodes.size() > 1) {
      std::pop_heap(nodes.begin(), nodes.end(), compare);
      auto node1 = nodes.back();
      nodes.pop_back();

      std::pop_heap(nodes.begin(), nodes.end(), compare);
      auto node2 = nodes.back();
      nodes.pop_back();

      nodes_owner.emplace_back(node1, node2);
      auto new_node = &nodes_owner.back();
      std::push_heap(nodes.begin(), nodes.end(), compare);
      nodes.push_back(new_node);
    }
    nodes.front()->dfs(lengths);
  }

  return HuffmanTree(lengths);
}


std::tuple<size_t,size_t, size_t> calculate_length_code(size_t length) {
  for (int i = 0; i < sizeof(length_table)/sizeof(length_table[0]); i++) {
    if (length_table[i][1] <= length && length < length_table[i+1][1]) {
      return {i + 257, length_table[i][0], length - length_table[i][1]};
    }
  }
  assert(0);
}

// dist_code, length of extra bits, content of extra bits
std::tuple<size_t, size_t, size_t> calculate_dist_code(size_t dist) {
  for (int i = 0; i < sizeof(distance_table)/sizeof(distance_table[0]); i++) {
    if (distance_table[i][1] <= dist && dist < distance_table[i+1][1]) {
      return {i, distance_table[i][0], dist - distance_table[i][1]};
    }
  }
  assert(0);
}

class Deflate {
 public:
  explicit Deflate(std::istream &is, size_t window_size) : is_(is), window_(window_size) { }

  void write_dynamic_huffman_header(
      HuffmanTree &lit_tree,
      HuffmanTree &dist_tree,
      BitStreamWriter &writer) {

    // HLIT
    writer.write_le(lit_tree.size()-257, 5);
    // HDIST
    writer.write_le(dist_tree.size()-1, 5);
    // HCLEN
    writer.write_le(19-4, 4);

    // No compression: all code lengths are 5 bits(use lit value)
    for (int i = 0; i < 19; i++) {
      writer.write_be(5, 3);
    }

    for (int i = 0; i < lit_tree.size(); i++) {
      auto [bits, nbits] = lit_tree.at(i);
      writer.write_be(nbits, 5);
    }

    for (int i = 0; i < dist_tree.size(); i++) {
      auto [bits, nbits] = dist_tree.at(i);
      writer.write_be(nbits, 5);
    }
  }

  // returns original size
  size_t deflate(std::ostream &os) {
    BitStreamWriter writer(os);
    size_t total_size = 0;
    while (is_.peek() != EOF) {
      lz77_encoded.clear();
      auto block_size = lz77_encode_block(65536);
      total_size += block_size;

      /*
       * Calculate count of lit/length and distance
       */
      std::vector<size_t> lit_counts(287, 0), dist_counts(30, 0);
      for (auto [c, distance, length] : lz77_encoded) {
        if (length != 0) {
          auto [lit_code, _1, _2] = calculate_length_code(length);
          auto [dist_code, _3, _4] = calculate_dist_code(distance);
          lit_counts[lit_code]++;
          dist_counts[dist_code]++;
        }
        if (c != EOF) {
          lit_counts[c]++;
        } else {
          lit_counts[256]++;
        }
      }

      // Create two Huffman trees
      auto lit_tree = create_huffman_tree(lit_counts);
      auto dist_tree = create_huffman_tree(dist_counts);

      // bfinal = false
      writer.write_le(0, 1);

      // 0b10: dynamic Huffman tree
      writer.write_le(0b10, 2);

      write_dynamic_huffman_header(lit_tree, dist_tree, writer);

      // the actual data
      for (auto [c, distance, length] : lz77_encoded) {
        if (length != 0) {
          {
            auto [len_code, nbits, bits] = calculate_length_code(length);
            lit_tree.write_value(writer, len_code);
          }

          {
            auto [dist_code, nbits, bits] = calculate_dist_code(distance);
            dist_tree.write_value(writer, dist_code);
            writer.write_be(bits, nbits);
          }
        }
        if (c == EOF) {
          lit_tree.write_value(writer, 256);
        } else {
          lit_tree.write_value(writer, c);
        }
      }
    }

    // This is the final block.
    // We use a non-compression empty block to indicate EOF
    // bfinal = true
    writer.write_le(1, 1);
    // no compression
    writer.write_le(0, 2);

    writer.flush_all();

    // empty data
    writer.write_le(0, 16);
    writer.write_le(0xffff, 16);

    return total_size;
  }

  // longest substring in the sliding window that matches current prefix
  // This function also manipulates the sliding window.
  // returns offset and length
  std::tuple<int, size_t, size_t> longest_prefix() {
    std::vector<size_t> matches, new_matches;
    bool first = true;
    size_t length = 0;
    while (true) {
      auto c = is_.get();

      if (first) {
        if (c == EOF) {
          return {c, 0, 0};
        }
        // offset in new_matches is the next offset to check
        window_.iterate([&new_matches, c = (uint8_t)c](uint8_t prev_c, size_t offset) {
          if (c == prev_c) {
            new_matches.push_back(offset);
          }
          return true;
        });
        first = false;
        if (new_matches.empty()) {
          window_.push_back(c);
          return std::make_tuple(c, 0, 0);
        }
      } else {
        if (c != EOF) {
          for (auto offset : matches) {
            if (window_.at(offset) == c) {
              new_matches.push_back(offset);
            }
          }
        }
        // including c == EOF
        if (new_matches.empty()) {
          assert(!matches.empty());
          // always use nearest match
          auto offset = matches.back();
          auto distance = window_.size() - offset;
          window_.push_back(c);
          return std::make_tuple(c, distance, length);
        }
      }

      auto popped = window_.push_back(c);
      if (!popped.has_value()) {
        // when no value is popped while pushing the value 'c', the new match index should increase 1
        for (auto &m : new_matches) {
          m++;
        }
      }

      std::swap(new_matches, matches);
      new_matches.clear();
      length++;
    }
  }

  // lz77 encode
  // yields lit, length, distance
  // reads from is_ and write to writer_
  size_t lz77_encode_block(size_t max_size) {
    constexpr int min_prefix_length = 3;
    constexpr int max_prefix_length = 4095;
    size_t size = 0;

//    std::cerr << "lz77 encode block" << std::endl;
    while (lz77_encoded.size() <= max_size) {

      // get longest prefix
      auto [c, distance, length] = longest_prefix();

      if (min_prefix_length <= length && length <= max_prefix_length) {
        // discard length from input
        lz77_encoded.emplace_back(c, distance, length);
//        std::cerr << "match(" << distance << "," << length << ")";
//        std::cout << (char)c;
      } else {
        // discard length from input
        if (length > 0) {
          for (size_t i = 0; i < length; i++) {
            auto ch = window_.at(window_.size() - length + i - 1);
            lz77_encoded.emplace_back(ch, 0, 0);
//            std::cout << (char)ch;
          }
        }
//        std::cout << (char)c;
        lz77_encoded.emplace_back(c, 0, 0);
      }
      size += length;
      if (c != EOF) {
        size += 1;
      } else {
        break;
      }
    }
    return size;
  }

 private:
  RingBuffer<uint8_t> window_;
  std::list<std::tuple<int, size_t, size_t>> lz77_encoded;
  std::istream &is_;
};

void compress_stream(std::istream &is, std::ostream &os) {
  // ID1
  os.put(0x1f);
  // ID2
  os.put(0x8b);
  // CM: deflate
  os.put(0x08);

  // flags:
  os.put(0);

  // uint32_t mtime
  os.put(0);
  os.put(0);
  os.put(0);
  os.put(0);

  // extra flag
  os.put(0);

  // OS: Linux
  os.put(0x03);

  Deflate deflate(is, 4096);
  auto original_size = deflate.deflate(os);

  // TODO: calculate CRC32
  uint32_t crc32 = 0;
  write_uint32(os, crc32);

  write_uint32(os, original_size);
}


int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "invalid arguements" << std::endl;
    return 1;
  }
  std::ostream *os = nullptr;
  std::istream *is = nullptr;

  std::string option = argv[1];
  std::string input_file_path = argv[2];
  std::string output_file_path = argv[3];
  std::ifstream ifs;
  std::ofstream ofs;

  if (input_file_path == "-") {
    is = &std::cin;
  } else {
    is = &ifs;
    ifs.open(input_file_path);
    if (!ifs) {
      std::cerr << "cannot open file '" << input_file_path << "'" << std::endl;
      return 1;
    }
  }

  if (output_file_path == "-") {
    os = &std::cout;
  } else {
    os = &ofs;
    ofs.open(output_file_path);
    if (!ofs) {
      std::cerr << "cannot open output file '" << output_file_path << "'" << std::endl;
      return 1;
    }
  }

  if (option.starts_with('d')) {
    decompress_stream(*is, *os);
  } else if (option.starts_with('c')) {
    compress_stream(*is, *os);
  } else {
    std::cerr << "Unknown options: '" << option << "'" << std::endl;
    return 1;
  }

  return 0;
}