#include <cassert>

#include <tuple>
#include <iomanip>
#include <iostream>
#include <fstream>

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

size_t length_table[][2] = {
    {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9}, {0, 10}, {1, 11}, {1, 13}, {1, 15}, {1, 17}, {2, 19}, {2, 23}, {2, 27}, {2, 31}, {3, 35}, {3, 43}, {3, 51}, {3, 59}, {4, 67}, {4, 83}, {4, 99}, {4, 115}, {5, 131}, {5, 163}, {5, 195}, {5, 227}, {0, 258}
};

// code -> (extra_bits, start_offset)
size_t distance_table[][2] = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 5}, {1, 7}, {2, 9}, {2, 13}, {3, 17}, {3, 25}, {4, 33}, {4, 49}, {5, 65}, {5, 97}, {6, 129}, {6, 193}, {7, 257}, {7, 385}, {8, 513}, {8, 769}, {9, 1025}, {9, 1537}, {10, 2049}, {10, 3073}, {11, 4097}, {11, 6145}, {12, 8193}, {12, 12289}, {13, 16385}, {13, 24577},{0,0}, {0,0}
};

class Deflate {
 public:
  explicit Deflate(std::istream& is) :is_(is), cached_byte_size(0), cached_byte(0) { }

  size_t deflate(std::ostream &os) {
    size_t total_size = 0;
    while (true) {
      auto [size, eof] = deflate_block(os);
      total_size += size;
      std::cout << "block size " << size<< std::endl;
      if (eof) {
        break;
      }
    }
    return total_size;
  }

  std::tuple<size_t, bool> deflate_block(std::ostream &os) {
    auto bfinal = read_bit(1);

    // 32KiB maximum
    auto block_buffer = new uint8_t[33 * 1024];
    size_t block_offset = 0;
    // assume static huffman
    auto btype = read_bit(2);
    assert(btype == 0b10);

    size_t original_size = 0;
    while (true) {
      auto code = read_fixed();
      assert(code <= 287);
      if (code < 256) {
        block_buffer[block_offset++] = code;
      } else if (code == 256) {
        break;
      } else {
        // first read length
        size_t length = 0, distance = 0;
        auto length_code = code - 257;
        {
          auto [extra_bits, start_off] = length_table[length_code];
          length = read_bit(extra_bits) + start_off;
        }

        // then distance
        auto dist_code = read_bit(5);
        {
          auto [ebits, start_off] = distance_table[dist_code];
          distance = read_bit(ebits) + start_off;
        }

        if (distance > block_offset) {
          throw InvalidFormat("distance should not > block_offset");
        }
        if (distance == 0) {
          throw InvalidFormat("distance should not be 0");
        }
        for (size_t i = 0; i < length; i++) {
          block_buffer[block_offset + i] = block_buffer[block_offset - distance + i];
        }
        block_offset += length;
      }
    }

    auto block_size = block_offset;

    original_size += block_size;
    os.write((char*)block_buffer, block_size);
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

  uint32_t read_fixed() {
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
  auto original_size = read_uint32(is);
//  if (decompressed_size != original_size) {
//    throw InvalidFormat("data corruption, original_size != decompressed_size");
//  }
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