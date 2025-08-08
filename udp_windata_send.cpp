// udp_windata_send.cpp
// C++ translation of OBS1.py (originally in Python)
// Sends one data frame per second via UDP.
// Data format follows custom "win" binary structure.
// Output is selectable between:
//   (1) UDP transmission to specified IP/port (default: 127.0.0.1:8001)
//   (2) Optional local file output (binary/text)
// Channel list is defined statically below:

#include <chrono>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#ifdef _WIN32
#include <winsock2.h>
#pragma comment(lib, "ws2_32.lib")
#else
#include <arpa/inet.h>
#include <sys/socket.h>
#include <unistd.h>
#endif

const char *UDP_IP = "127.0.0.1";
const int UDP_PORT = 8001;
const char *SAVE_RAW_FILE = "data.win";
const char *SAVE_HEX_FILE = "data.txt";

// ======== CHANNEL LIST CONFIGURATION ==========
const std::vector<uint16_t> channels = {
    0xd800, 0xd801, 0xd802, 0xd805 // Add or remove channels here
};

uint8_t bcd_encode(int value) { return ((value / 10) << 4) | (value % 10); }

std::vector<uint8_t> get_bcd_time_bytes() {
    std::time_t t = std::time(nullptr);
    std::tm *now = std::localtime(&t);
    return {bcd_encode(now->tm_year % 100), bcd_encode(now->tm_mon + 1), bcd_encode(now->tm_mday),
            bcd_encode(now->tm_hour),       bcd_encode(now->tm_min),     bcd_encode(now->tm_sec)};
}

std::vector<uint8_t> build_frame_second(uint8_t packet_no, uint8_t code,
                                        const std::vector<uint16_t> &channels, int sample_size = 1,
                                        int sample_rate = 2) {
    std::vector<uint8_t> frame;
    frame.push_back(packet_no);
    frame.push_back(packet_no);
    frame.push_back(code);

    int blosize = 5 + 6 + channels.size() * (4 + 4 + sample_rate - 1);
    uint16_t blosize_short = static_cast<uint16_t>(blosize - 3);
    frame.push_back(blosize_short >> 8);
    frame.push_back(blosize_short & 0xFF);

    auto bcd_time = get_bcd_time_bytes();
    frame.insert(frame.end(), bcd_time.begin(), bcd_time.end());

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 5);

    for (uint16_t ch : channels) {
        frame.push_back(ch >> 8);
        frame.push_back(ch & 0xFF);
        uint16_t sr_field = ((sample_size & 0xF) << 12) | (sample_rate & 0xFFF);
        frame.push_back(sr_field >> 8);
        frame.push_back(sr_field & 0xFF);
        uint32_t first_sample = dis(gen);
        frame.push_back((first_sample >> 24) & 0xFF);
        frame.push_back((first_sample >> 16) & 0xFF);
        frame.push_back((first_sample >> 8) & 0xFF);
        frame.push_back(first_sample & 0xFF);
        uint8_t diff_sample = dis(gen);
        frame.push_back(diff_sample);
    }

    return frame;
}

std::string hexdump(const std::vector<uint8_t> &data) {
    std::ostringstream oss;
    for (size_t i = 0; i < data.size(); i += 16) {
        oss << std::setw(8) << std::setfill('0') << std::hex << static_cast<int>(i) << "  ";
        for (size_t j = 0; j < 16 && i + j < data.size(); j += 2) {
            if (i + j + 1 < data.size()) {
                uint16_t word = (data[i + j] << 8) | data[i + j + 1];
                oss << std::setw(4) << std::setfill('0') << std::hex << word << " ";
            } else {
                oss << std::setw(2) << std::setfill('0') << std::hex
                    << static_cast<int>(data[i + j]) << "  ";
            }
        }
        oss << "\n";
    }
    return oss.str();
}

int main() {
#ifdef _WIN32
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
#endif

    std::string save_binary_choice, save_text_choice;
    std::cout << "Save binary output to file ? (y/n): ";
    std::cin >> save_binary_choice;
    std::cout << "Save hexdump output to file ? (y/n): ";
    std::cin >> save_text_choice;

    bool save_binary = (save_binary_choice == "y" || save_binary_choice == "Y");
    bool save_text = (save_text_choice == "y" || save_text_choice == "Y");

    int sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(UDP_PORT);
    addr.sin_addr.s_addr = inet_addr(UDP_IP);

    uint8_t packet_no = 0;
    while (true) {
        auto frame = build_frame_second(packet_no, 0xA0, channels);
        sendto(sockfd, reinterpret_cast<const char *>(frame.data()), frame.size(), 0,
               (sockaddr *)&addr, sizeof(addr));
        std::cout << "Sent frame, size = " << frame.size() << " bytes\n";

        if (save_binary) {
            std::ofstream raw_out(SAVE_RAW_FILE, std::ios::app | std::ios::binary);
            raw_out.write(reinterpret_cast<const char *>(frame.data()), frame.size());
        }

        if (save_text) {
            std::ofstream hex_out(SAVE_HEX_FILE, std::ios::app);
            hex_out << "Packet#" << static_cast<int>(packet_no) << ", size=" << frame.size()
                    << "\n";
            hex_out << hexdump(frame) << "\n";
        }

        std::this_thread::sleep_for(std::chrono::seconds(1));
        packet_no = (packet_no + 1) % 256;
    }

#ifdef _WIN32
    closesocket(sockfd);
    WSACleanup();
#else
    close(sockfd);
#endif
    return 0;
}
