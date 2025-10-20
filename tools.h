
#include <sstream>

std::string join(const std::vector<std::string>& vec) {
        std::string tmp;
        for (size_t i = 0; i < vec.size(); ++i) {
                tmp  +=vec[i];
        }
        tmp+='@';
        return tmp;
}



std::vector<std::string> split(const std::string& str) {
        std::vector<std::string> result;
        std::string temp;
        for (char ch : str) {
                if (ch == '@')
                {
                        result.push_back(temp);
                        temp.clear();
                }
                else
                {
                        temp += ch;
                }
        }

        result.push_back(temp);
        return result;
}


std::string read_file(std::string fileName) {

std::ifstream file(fileName);

if (!file) {
	std::cerr << "Can't open file." << std::endl;
	return "1";
}

	std::stringstream buffer;
	buffer << file.rdbuf();
	return  buffer.str();

}
