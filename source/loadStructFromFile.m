function result = loadStructFromFile(fileName, environmentName)
tmp = load(fileName, environmentName);
result = tmp.(environmentName);