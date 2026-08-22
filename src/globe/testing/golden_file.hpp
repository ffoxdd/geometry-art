#ifndef GLOBEART_SRC_GLOBE_TESTING_GOLDEN_FILE_HPP_
#define GLOBEART_SRC_GLOBE_TESTING_GOLDEN_FILE_HPP_

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <string>

namespace globe::testing {

// Compares generated text against a file committed beside the test. Set
// UPDATE_GOLDEN=1 to rewrite them after an intended change, so the diff of
// what the pipeline emits shows up in review rather than being asserted
// piecemeal.
class GoldenFile {
 public:
    explicit GoldenFile(std::string relative_path);

    void expect_matches(const std::string& content) const;

 private:
    std::string _path;

    [[nodiscard]] static bool updating();
    void rewrite(const std::string& content) const;
    [[nodiscard]] std::string read() const;
};

inline GoldenFile::GoldenFile(std::string relative_path) :
    _path(std::string(GLOBE_SOURCE_DIR) + "/" + std::move(relative_path)) {
}

inline void GoldenFile::expect_matches(const std::string& content) const {
    if (updating()) {
        rewrite(content);
        GTEST_SKIP() << "rewrote " << _path;
    }

    ASSERT_TRUE(std::filesystem::exists(_path))
        << _path << " is missing; run with UPDATE_GOLDEN=1 to create it";

    EXPECT_EQ(content, read())
        << "output differs from " << _path << "; run with UPDATE_GOLDEN=1 if intended";
}

inline bool GoldenFile::updating() {
    const char* flag = std::getenv("UPDATE_GOLDEN");
    return flag != nullptr && flag[0] != '\0';
}

inline void GoldenFile::rewrite(const std::string& content) const {
    std::filesystem::create_directories(std::filesystem::path(_path).parent_path());
    std::ofstream stream(_path);
    stream << content;
}

inline std::string GoldenFile::read() const {
    std::ifstream stream(_path);
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}

} // namespace globe::testing

#endif //GLOBEART_SRC_GLOBE_TESTING_GOLDEN_FILE_HPP_
