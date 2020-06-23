#pragma once
#include <string>
#include <map>
#include <JuceHeader.h>


namespace JuceUTF16toStdString
{

	static const std::map<unsigned int, char> my_map = {
	{ 225, '�' },
	{ 269, '�' },
	{ 271, '�' },
	{ 233, '�' },
	{ 283, '�' },
	{ 237, '�' },
	{ 328, '�' },
	{ 243, '�' },
	{ 345, '�' },
	{ 353, '�' },
	{ 357, '�' },
	{ 250, '�' },
	{ 367, '�' },
	{ 253, '�' },
	{ 382, '�' },
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
