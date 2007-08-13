def colorize(text, color_code)
  # Make sure the terminal supports color
  if (ENV.has_key?("COLORTERM") and ENV["COLORTERM"] == 1) or
    (ENV.has_key?("TERM") and ENV["TERM"].include?("color"))
    "#{color_code}#{text}\e[0m"
  else
    text
  end
end

def red(text);     colorize(text, "\e[31m"); end
def green(text);   colorize(text, "\e[32m"); end
def blue(text);    colorize(text, "\e[34m"); end
def yellow(text);  colorize(text, "\e[33m"); end
def magenta(text); colorize(text, "\e[35m"); end
def cyan(text);    colorize(text, "\e[36m"); end
def white(text);   colorize(text, "\e[37m"); end
def bold(text);    colorize(text, "\e[1m");  end
