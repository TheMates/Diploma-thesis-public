#pragma once
#include <string>
#include <map>
#include <JuceHeader.h>


namespace JuceUTF16toStdString
{

	static const std::map<unsigned int, char> my_map = {
	{ 225, 'á' },
	{ 269, 'è' },
	{ 271, 'ï' },
	{ 233, 'é' },
	{ 283, 'ì' },
	{ 237, 'í' },
	{ 328, 'ò' },
	{ 243, 'ó' },
	{ 345, 'ø' },
	{ 353, 'š' },
	{ 357, '' },
	{ 250, 'ú' },
	{ 367, 'ù' },
	{ 253, 'ı' },
	{ 382, '' },
	};

	class Converter
	{
	public:
		static std::string Convert(juce::CharPointer_UTF16& utf16String)
		{
			int size = utf16String.length();
			std::string outString = "";
			char ch;
			for (auto i = 0; i < size; ++i)
			{
				auto code = utf16String[i];
				if (code > 127)
					ch = my_map.at(code);
				else
					ch = code;
				outString += ch;
			}
			return outString;
		}
	};
}
