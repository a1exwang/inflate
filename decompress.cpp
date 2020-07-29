#include <cassert>

#include <tuple>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

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

constexpr size_t length_table[][2] = {
    {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9}, {0, 10}, {1, 11}, {1, 13}, {1, 15}, {1, 17}, {2, 19}, {2, 23}, {2, 27}, {2, 31}, {3, 35}, {3, 43}, {3, 51}, {3, 59}, {4, 67}, {4, 83}, {4, 99}, {4, 115}, {5, 131}, {5, 163}, {5, 195}, {5, 227}, {0, 258}
};

// code -> (extra_bits, start_offset)
constexpr size_t distance_table[][2] = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 5}, {1, 7}, {2, 9}, {2, 13}, {3, 17}, {3, 25}, {4, 33}, {4, 49}, {5, 65}, {5, 97}, {6, 129}, {6, 193}, {7, 257}, {7, 385}, {8, 513}, {8, 769}, {9, 1025}, {9, 1537}, {10, 2049}, {10, 3073}, {11, 4097}, {11, 6145}, {12, 8193}, {12, 12289}, {13, 16385}, {13, 24577},{0,0}, {0,0}
};

class Deflate {
 public:
  explicit Deflate(std::istream& is) :is_(is), cached_byte_size(0), cached_byte(0) { }

  size_t deflate(std::ostream &os) {
    size_t total_size = 0;
    while (true) {
      std::cout << "start to deflate block off=" << total_size << std::endl;
      auto [size, eof] = deflate_block(os);
      total_size += size;
      std::cout << "block size " << size<< std::endl;
      if (eof) {
        break;
      }
    }
    return total_size;
  }

  struct HuffmanTree {
    HuffmanTree() { }
    explicit HuffmanTree(const std::vector<size_t> &code_lengths)
        :counts(code_lengths.size(), 0), offsets(code_lengths.size(), 0), skip(code_lengths.size(), 0) {

      std::vector<std::vector<int>> n_codes(code_lengths.size(), std::vector<int>(0));
      for (int code = 0; code < code_lengths.size(); code++) {
        counts[code_lengths[code]]++;
        n_codes[code_lengths[code]].push_back(code);
      }

      counts[0] = 0;
      for (int i = 1; i < offsets.size(); i++) {
        offsets[i] = counts[i-1] + offsets[i-1];
      }

      uint64_t last_skip = 0;
      for (int i = 1; i < counts.size(); i++) {
        last_skip = (last_skip + counts[i-1]) << 1u;
        skip[i] = last_skip;
      }

      for (int i = 1; i < n_codes.size(); i++) {
        for (auto code : n_codes[i]) {
          symbols.push_back(code);
        }
      }
    }
    void print() {
      for (int i = 0; i < symbols.size(); i++) {
        std::cout << "symbols[" << i << "] = " << symbols[i] << std::endl;
      }
      for (auto offset : offsets) {
        std::cout << offset << " ";
      }
      std::cout << std::endl;
      for (auto count : counts) {
        std::cout << count << " ";
      }
      std::cout << std::endl;

    }
    std::vector<size_t> offsets;
    std::vector<size_t> counts;
    std::vector<size_t> symbols;
    std::vector<uint64_t> skip;
  };

  uint32_t read_value_from_huffman(const HuffmanTree &tree) {
    uint64_t value = 0;
    size_t bit_length = 0;
    while (true) {
      value |= read_bit(1);
      bit_length++;
      assert(bit_length < tree.offsets.size());

      auto off = value - tree.skip[bit_length];
      if (off < tree.counts[bit_length]) {
        auto index = tree.offsets[bit_length] + off;
        return tree.symbols[index];
      }
      value <<= 1u;
    }
  }

  std::vector<size_t> generate_huffman_table(HuffmanTree &huffman_tree, const size_t lit_count) {
    std::vector<size_t> huffman_table;
    int prev = -1;

    for (int lit = 0; lit < lit_count; ) {
      auto value = read_value_from_huffman(huffman_tree);
//      std::cout << "read " << (int)value << std::endl;

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
        assert(prev != -1);
        int copy_previous_count = read_bit_le(2) + 3;
        for (int i = 0; i < copy_previous_count; i++) {
          huffman_table.push_back(prev);
        }
        lit += copy_previous_count;
      } else if (17 <= value && value <= 18) {
        int zero_count = 0;
        if (value == 17) {
          zero_count = read_bit_le(3) + 3;
        } else {
          zero_count = read_bit_le(7) + 11;
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

  std::tuple<HuffmanTree, HuffmanTree> read_dynamic_huffman_tree_header() {
    const auto hlit = read_bit_le(5) + 257;
    const auto hdist = read_bit_le(5) + 1;
    const auto hclen = read_bit_le(4) + 4;

    std::vector<size_t> code_lengths(19, 0);

    int index[] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
    for (int i = 0; i < hclen; i++) {
      auto code_length = read_bit_le(3);
      code_lengths[index[i]] = code_length;
    }

    HuffmanTree code_length_huffman_tree(code_lengths);
//    code_length_huffman_tree.print();

    auto lit_table = generate_huffman_table(code_length_huffman_tree, hlit);
    HuffmanTree lit_tree(lit_table);

    auto dist_table = generate_huffman_table(code_length_huffman_tree, hdist);
    HuffmanTree dist_tree(dist_table);
//    dist_huffman_tree->print();
    return std::make_tuple(lit_tree, dist_tree);
  }

  size_t deflate_huffman(
      std::ostream &os,
      bool fixed,
      const HuffmanTree &lit_tree,
      const HuffmanTree &dist_tree) {

    size_t block_offset = 0;
    size_t original_size = 0;
    while (true) {

      uint32_t code;
      if (fixed) {
        code = read_fixed_lit();
      } else {
        code = read_value_from_huffman(lit_tree);
      }
      assert(code <= 287);
      if (code < 256) {
        block_buffer[block_offset++] = code;
//        std::cerr <<(char)code;
//        std::cout << "decode lit 0x" << std::hex << code << std::endl;
      } else if (code == 256) {
//        std::cout << std::endl << "---------------------------" << std::endl;
        break;
      } else {
        // first read length
        size_t length = 0, distance = 0;
        auto length_code = code - 257;
        {
          auto [extra_bits, start_off] = length_table[length_code];
          length = read_bit_le(extra_bits) + start_off;
        }

        uint32_t dist_code = 0;
        if (fixed) {
          // then distance
          dist_code = read_bit(5);
        } else {
          dist_code = read_value_from_huffman(dist_tree);
//          std::cout << "huff " << dist_huffman_tree->symbols[0] <<  " dist_code: " << dist_code << std::endl;
        }

        {
          auto [extra_bits, start_off] = distance_table[dist_code];
          distance = read_bit_le(extra_bits) + start_off;
        }

        if (distance == 0) {
          throw InvalidFormat("distance should not be 0");
        }
        if (distance > block_offset) {
          // FIXME: Currently we only support "distance" within one block distance
          assert(prev_block_size > (distance - block_offset));
          auto prev_block_off = prev_block_size - (distance - block_offset);
          if (length <= prev_block_size - prev_block_off) {
            std::copy(&prev_block_buffer[prev_block_off], prev_block_buffer.data()+prev_block_off+length, &block_buffer[block_offset]);
          } else {
            std::copy(&prev_block_buffer[prev_block_off], prev_block_buffer.data()+prev_block_size, &block_buffer[block_offset]);

            auto length_copied = (prev_block_size - prev_block_off);
            auto remaining_length = length - length_copied;
            for (int i = 0; i < remaining_length; i++) {
              block_buffer[block_offset + length_copied + i] = block_buffer[i];
            }
          }
        } else {
          // length might <= distance according to the standard.
          // The correct behavior is as follows
          assert(block_offset + length <= block_buffer.size());
          for (int i = 0; i < length; i++) {
            block_buffer[block_offset + i] = block_buffer[block_offset - distance + i];
          }
        }
//        std::copy(start, start + length, block_buffer.data() + block_offset);

//        std::cout << "match " << length << ", " << distance;
//        std::cerr.write((char*)block_buffer.data()+block_offset, length);
//        std::cerr << std::endl;

        block_offset += length;
      }
    }

    auto block_size = block_offset;

    original_size += block_size;
    os.write((char*)block_buffer.data(), block_size);
    return block_size;
  }


  template <typename T>
  T read_le() {
    T value;
    for (int i = 0; i < sizeof(T); i++) {
      reinterpret_cast<uint8_t*>(&value)[i] = is_.get();
      assert(is_);
    }
    return value;
  }

  std::tuple<size_t, bool> deflate_block(std::ostream &os) {
    auto bfinal = read_bit(1);

    // 0b10: static huffman
    // 0b01: dynamic huffman
    auto btype = read_bit(2);
    size_t block_size = 0;
    if (btype == 0b10) {
      block_size = deflate_huffman(os, true, HuffmanTree(), HuffmanTree());
    } else if (btype == 0b01) {
      auto [lit_tree, dist_tree] = read_dynamic_huffman_tree_header();
      block_size = deflate_huffman(os, false, lit_tree, dist_tree);
    } else if (btype == 0) {
      discard_remaining_bits();
      assert(cached_byte_size == 0);

      auto len = read_le<uint16_t>();
      auto nlen = read_le<uint16_t>();
      assert((len ^ nlen) == 0xffff);
      for (int i = 0; i < len; i++) {
        os.put(is_.get());
        assert(is_);
      }

      block_size = len;
    } else {
      assert(0);
    }

    std::swap(block_buffer, prev_block_buffer);
    prev_block_size = block_size;
    return {block_size, bfinal};
  }

 private:
  void ensure_cached_byte() {
    if (cached_byte_size == 0) {
      cached_byte_size = 8;
      auto c = is_.get();
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
    assert(n < 32);
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
    assert(n < 32);
    uint32_t result = 0;
    for (size_t i = 0; i < n; i++) {
      ensure_cached_byte();
      auto target_bit = ((uint32_t)cached_byte >> (8-cached_byte_size)) & 1u;
      result |= (target_bit << i);
      cached_byte_size--;
    }
    return result;
  }

  uint32_t read_fixed_lit() {
    auto prefix4 = read_bit(4);
    if (prefix4 < 0b0011) {
      // 256 - 279
      auto part3 = read_bit(3);
      if (prefix4 & 1u) {
        part3 += 8;
      }
      return part3 + 256;
    } else if (prefix4 >= 0b0011 && prefix4 < 0b1100) {
      auto prefix5 = (prefix4 << 1u) | read_bit(1);
      if (prefix5 >= 0b00110 && prefix5 < 0b11000) {
        // 0 - 143
        return ((prefix5 << 3u) | read_bit(3)) - 0x30;
      } else if (prefix5 >= 0b11000 && prefix5 < 0b11001) {
        // 280 - 287
        return read_bit(3) + 280;
      } else {
        // 144 - 255, code starts from 0b110010000
        return ((prefix5 << 4u) | read_bit(4)) - 0b110010000 + 144;
      }
    } else {
      assert(0);
    }
  }

 private:
  // 32KiB maximum
  std::vector<uint8_t> block_buffer = std::vector<uint8_t>(1024*1024);
  std::vector<uint8_t> prev_block_buffer = std::vector<uint8_t>(1024*1024);
  size_t prev_block_size = 0;

  uint8_t cached_byte;
  size_t cached_byte_size;
  std::istream& is_;
};

// RFC1952
void read_stream(std::istream &is, std::ostream &os) {

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
  // not implemented
  assert(flags == 0);

  // mtime
  uint32_t mtime = read_uint32(is);
  {
    std::time_t temp = mtime;
    std::tm* t = std::gmtime(&temp);
    std::stringstream ss; // or if you're going to print, just input directly into the output stream
    ss << std::put_time(t, "%Y-%m-%d %I:%M:%S %p");
    std::string mtime_string = ss.str();
    std::cout << "mtime: " << mtime_string << std::endl;
  }

  // extra flag
  skip(is, 1);

  // OS
  skip(is, 1);

  Deflate deflate(is);
  size_t decompressed_size = deflate.deflate(os);

  auto crc32 = read_uint32(is);
  std::cout << "crc32: " << crc32 << std::endl;
  auto original_size = read_uint32(is);
  if (decompressed_size != original_size) {
    throw InvalidFormat("data corruption, original_size != decompressed_size");
  }
  std::cout << "original_size: " << original_size << std::endl;
}


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "invalid arguements" << std::endl;
    return 1;
  }

  std::string input_file_path = argv[1];
  std::string output_file_path = argv[2];
  std::ifstream ifs(input_file_path);
  if (!ifs) {
    std::cerr << "cannot open file '" << input_file_path << "'" << std::endl;
    return 1;
  }
  std::ofstream ofs(output_file_path);
  if (!ofs) {
    std::cerr << "cannot open output file '" << output_file_path << "'" << std::endl;
    return 1;
  }
  read_stream(ifs, ofs);
  return 0;
}