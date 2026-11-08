#ifndef RK4_IO_CONFIG_H
#define RK4_IO_CONFIG_H

//----------------------------------------------------------------------
// 极简配置解析（key = value，支持 [section]）
//
// 不引入 JSON/YAML 库：格式简单、可读、够用。
//
//   [reaction]
//   beam_A = 12
//   beam_Z = 4
//
//   [beam]
//   ke_mev = 153.9
//
// 读取时用 "reaction.beam_A" 这样的键。
//----------------------------------------------------------------------

#include <map>
#include <string>

namespace rk4 {

class Config {
public:
    bool Load(const std::string& path);

    bool GetBool(const std::string& key,bool def = false) const;
    bool Has(const std::string& key) const;
    std::string GetString(const std::string& key, const std::string& def = "") const;
    double GetDouble(const std::string& key, double def = 0.0) const;
    int GetInt(const std::string& key, int def = 0) const;

    const std::map<std::string, std::string>& Values() const { return values_; }

private:
    std::map<std::string, std::string> values_;
};

} // namespace rk4

#endif // RK4_IO_CONFIG_H
