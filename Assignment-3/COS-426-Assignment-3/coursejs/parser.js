// this file allows us to use commands in the url for batch processing
// this is some crazy magic. don't read too closely.

var Parser = Parser || {};

// gets commands from URL and separates into an array of objects
Parser.getCommands = function(url) {
  // eliminate anything in string up through '?'
  var pos = url.lastIndexOf("?");
  url = url.substr(pos + 1);

  // split the URL by the word 'apply'
  url = url.replace(/apply/g, "|");
  var applies = url.split("|");
  var paramsArray = [];

  for (var i = 0; i < applies.length; i++) {
    var params = {};
    var commands = applies[i].split("&");
    for (var j = 0; j < commands.length; j++) {
      var command = commands[j];
      var parts = command.split("=");

      if (command.length == 0) {
        continue;
      } else if (parts.length == 1) {
        params[parts[0]] = true;
      } else {
        params[parts[0]] = parts[1];
      }
    }
    paramsArray.push(params);
  }
  return paramsArray;
};

// turns the commands from strings into objects
Parser.parseCommands = function(applies) {
  var result = [];
  for (var i = 0; i < applies.length; i++) {
    var commands = applies[i];
    var apply = {};
    // update values from url using parser
    for (var cmd in commands) {
      var v = unescape(commands[cmd]).replace("+", " ");
      if (!isNaN(parseFloat(v))) {
        v = parseFloat(v);
      } else if (v == "true") {
        v = true;
      } else if (v == "false") {
        v = false;
      }
      apply[cmd] = v;
    }
    result.push(apply);
  }
  return result;
};

Parser.parseNumbers = function(str) {
  var numbers = [];
  str.split(",").forEach(function(interval) {
    var parts = interval.split("-");
    if (parts.length == 1) {
      numbers.push(parseInt(parts[0]));
    } else {
      var start = parseInt(parts[0]);
      var end = parseInt(parts[1]);
      $.merge(numbers, Parser.range(start, end + 1));
    }
  });
  return numbers;
};

Parser.parseJson = function(jsonFile) {
  var request = new XMLHttpRequest();
  request.open("GET", jsonFile, false);
  request.overrideMimeType("application/json");
  request.send(null);
  var obj = JSON.parse(request.responseText);
  return obj;
};

Parser.parseTxt = function(textFile) {
  var request = new XMLHttpRequest();
  request.open("GET", textFile, false);
  request.overrideMimeType("text/plain");
  request.send(null);
  return request.responseText;
};
