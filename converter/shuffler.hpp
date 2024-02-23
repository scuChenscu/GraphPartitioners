#pragma once

#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <algorithm>
#include <random>
#include <utility>
#include "conversions.hpp"

class Shuffler : public Converter {
private:
    struct work_t {
        Shuffler *shuffler;
        int nchunks;
        std::vector<edge_t> chunk_buf;

        void operator()() {
            std::shuffle(chunk_buf.begin(), chunk_buf.end(),
                         std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
            int file =
                    open(shuffler->chunk_filename(nchunks).c_str(),
                         O_WRONLY | O_CREAT, S_IROTH | S_IWOTH | S_IWUSR | S_IRUSR);
            size_t chunk_size = chunk_buf.size() * sizeof(edge_t);
            writea(file, (char *) &chunk_buf[0], chunk_size);
            close(file);
            chunk_buf.clear();
        }
    };

    size_t chunk_bufsize;
    int nchunks;
    std::vector<edge_t> chunk_buf;

    std::string chunk_filename(int chunk);

    void chunk_clean();

    void cwrite(edge_t e, bool flush = false);

public:
    Shuffler(std::string input) : Converter(std::move(input)) {}

    bool done() { return is_exists(shuffled_binary_edgelist_name(filename)); }

    void init(int memorysize);

    void finalize();

    void add_edge(vid_t source, vid_t target);

    static void writea(int f, char *buf, size_t nbytes) {
        size_t nwritten = 0;
        while (nwritten < nbytes) {
            ssize_t a = write(f, buf, nbytes - nwritten);
            PCHECK(a != ssize_t(-1)) << "Could not write " << (nbytes - nwritten)
                                     << " bytes!";
            buf += a;
            nwritten += a;
        }
    }

};

