Here are the tables:
- studies: (id, datetime, name, study_params)
- datasets: (id, study_id, SimulationParams, data [stored on disk])
- runs: (id, dataset_id, model_name, results [stored on disk])

```mermaid
erDiagram
    direction LR
    __studies__ {
        int id PK
        str name
        datetime datetime
        StudyParams study_params
    }
    __sim_params__{
        int id PK
        int study_id FK
        SimulationParams simulation_params
    }
    __datasets__ {
        int id PK
        int sim_param_id FK
        COMBData data
    }
    __runs__ {
        int id PK
        int dataset_id FK
        str model_name
        InferenceData results
    }
    __studies__ one to zero or more __sim_params__ : contains
    __sim_params__ one to zero or more __datasets__ : contains
    __datasets__ one to zero or more __runs__ : contains
```