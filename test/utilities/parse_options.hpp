#include <cerrno>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

enum class ParseResult
{
  Success,
  Stop,
  Error
};

struct OptionSpec
{
  const char* shortFlag; // "-b"
  const char* longFlag;  // "--beta"
  bool takesValue;       // true if the option expects following values
  const char* help;      // help text (include defaults in text)
  std::function<ParseResult(const char* v)> set;
};

static void printUsage(const char* progname,
                       const std::vector<OptionSpec>& options)
{
  std::cout << "Usage: " << progname << " [options]\n\n";
  std::cout << "Options:\n";
  for (const auto& s : options)
  {
    std::cout << "  " << s.shortFlag << ", " << s.longFlag;
    if (s.takesValue) std::cout << " <value>";
    std::cout << "\n      " << (s.help ? s.help : "") << "\n";
  }
}

static int parseArguments(int argc, char* argv[], std::vector<OptionSpec>& options)
{
  // Add help flags so the example does not need to
  options.push_back(OptionSpec{"-h", "--help", false, "Print this help message",
                               [&](const char*)
                               {
                                 printUsage(argv[0], options);
                                 return ParseResult::Stop; // not an error
                               }});

  auto findSpec = [&](const char* arg) -> const OptionSpec*
  {
    for (const auto& s : options)
      if ((s.shortFlag && std::strcmp(arg, s.shortFlag) == 0) ||
          (s.longFlag && std::strcmp(arg, s.longFlag) == 0))
        return &s;
    return nullptr;
  };

  for (int i = 1; i < argc; ++i)
  {
    const char* arg        = argv[i];
    const OptionSpec* spec = findSpec(arg);
    if (!spec)
    {
      std::cerr << "Error: Unknown option '" << arg << "'\n\n";
      printUsage(argv[0], options);
      return 1;
    }

    const char* value = nullptr;
    if (spec->takesValue)
    {
      if (i + 1 >= argc)
      {
        std::cerr << "Error: " << arg << " requires a value\n\n";
        printUsage(argv[0], options);
        return 1;
      }
      value = argv[++i];
    }

    ParseResult r = spec->set(value);
    if (r == ParseResult::Stop) return 1;
    if (r == ParseResult::Error)
    {
      std::cerr << "Error: Invalid value for " << arg << "\n\n";
      printUsage(argv[0], options);
      return 1;
    }
  }

  return 0;
}

// Helpers
static inline ParseResult setBoolTrue(bool& out)
{
  out = true;
  return ParseResult::Success;
}

template<typename Real>
static ParseResult setReal(Real& out, const char* v)
{
  if (!v) return ParseResult::Error;
  out = static_cast<Real>(std::stold(v));
  return ParseResult::Success;
  //return true;
}

static inline ParseResult setString(std::string& out, const char* v)
{
  if (!v) return ParseResult::Error;
  out = v;
  return ParseResult::Success;
  //return true;
}
