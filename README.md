# modelica-vocabularies

## Build
Equivalent to owlReason task
```
./gradlew build
```

## Generate Docs
You must first have Bikeshed (the app itself) installed from [here](https://tabatkins.github.io/bikeshed/#install-final)
```
./gradlew generateDocs
```
Note: if bikeshed is not in the PATH, you can add -pBIKESHED=path/to/bikeshed argument

## Run OWL Reasoner
```
./gradlew owlReason
```

## Start Fuseki Server
```
./gradlew startFuseki
```

## Stop Fuseki Server
```
./gradlew stopFuseki
```

## Load to Fuseki Dataset
```
./gradlew owlLoad
```

## Run SPARQL Queries
```
./gradlew owlQuery
```

## Publish to Maven Local
```
./gradlew publishToMavenLocal
```
