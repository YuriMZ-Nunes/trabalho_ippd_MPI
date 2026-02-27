#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstddef>
#include <string>

const int W = 200;
const int H = 200;
const int T = 500;
const int S = 20;
const int N_AGENTS = 1000;

struct Cell {
    int tipo;
    int recurso;
    bool acessivel;
    int consumo_acumulado;
};

struct Agent {
    int id;
    int x, y;
    int energia;
};

struct Destino {
    bool eh_local;
    int processo_destino;
    int novo_x;
    int novo_y;
};

Destino decidir(Agent a, int offsetX, int local_W, int local_H, 
                const std::vector<std::vector<Cell>>& grid_local,
                const std::vector<Cell>& halo_esq,
                const std::vector<Cell>& halo_dir,
                int rank, int rank_esq, int rank_dir) {
    
    int lx = a.x - offsetX;
    int ly = a.y;

    int melhor_recurso = -1;
    int num_empatados = 0;
    Destino melhor_destino = {true, rank, a.x, a.y};

    unsigned int seed = a.id * 73856093 ^ a.x * 19349663 ^ a.y * 83492791;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            int nx = lx + dx;
            int ny = ly + dy;
            int global_x = a.x + dx;
            int global_y = a.y + dy;

            if (ny < 0 || ny >= local_H) continue;

            int rec = -1;
            bool acessivel = false;
            int dest_rank = rank;
            bool dest_local = true;

            if (nx < 0) {
                if (rank_esq != MPI_PROC_NULL) {
                    rec = halo_esq[ny].recurso;
                    acessivel = halo_esq[ny].acessivel;
                    dest_rank = rank_esq;
                    dest_local = false;
                }
            } else if (nx >= local_W) {
                if (rank_dir != MPI_PROC_NULL) {
                    rec = halo_dir[ny].recurso;
                    acessivel = halo_dir[ny].acessivel;
                    dest_rank = rank_dir;
                    dest_local = false;
                }
            } else {
                rec = grid_local[nx][ny].recurso;
                acessivel = grid_local[nx][ny].acessivel;
            }

            if (acessivel) {
                if (rec > melhor_recurso) {
                    melhor_recurso = rec;
                    num_empatados = 1;
                    melhor_destino.eh_local = dest_local;
                    melhor_destino.processo_destino = dest_rank;
                    melhor_destino.novo_x = global_x;
                    melhor_destino.novo_y = global_y;
                } 
                else if (rec == melhor_recurso) {
                    num_empatados++;
                    
                    if (rand_r(&seed) % num_empatados == 0) {
                        melhor_destino.eh_local = dest_local;
                        melhor_destino.processo_destino = dest_rank;
                        melhor_destino.novo_x = global_x;
                        melhor_destino.novo_y = global_y;
                    }
                }
            }
        }
    }
    return melhor_destino;
}

// --- Tipos Derivados MPI Globais ---
MPI_Datatype MPI_AGENT;
MPI_Datatype MPI_CELL;

void registrar_tipos_mpi() {
    int count_agent = 4;
    int blocklengths_agent[4] = {1, 1, 1, 1};
    MPI_Datatype types_agent[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint offsets_agent[4] = {offsetof(Agent, id), offsetof(Agent, x), offsetof(Agent, y), offsetof(Agent, energia)};
    MPI_Type_create_struct(count_agent, blocklengths_agent, offsets_agent, types_agent, &MPI_AGENT);
    MPI_Type_commit(&MPI_AGENT);

    int count_cell = 4;
    int blocklengths_cell[4] = {1, 1, 1, 1};
    MPI_Datatype types_cell[4] = {MPI_INT, MPI_INT, MPI_C_BOOL, MPI_INT};
    MPI_Aint offsets_cell[4] = {offsetof(Cell, tipo), offsetof(Cell, recurso), offsetof(Cell, acessivel), offsetof(Cell, consumo_acumulado)};
    MPI_Type_create_struct(count_cell, blocklengths_cell, offsets_cell, types_cell, &MPI_CELL);
    MPI_Type_commit(&MPI_CELL);
}

void executar_carga(int r) {
    int limit = r * 100; 
    volatile double dummy = 0.0;
    for (int c = 0; c < limit; ++c) {
        dummy += 1.0; 
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    registrar_tipos_mpi();

    int local_W = W / size; 
    int local_H = H;
    int offsetX = rank * local_W;
    int offsetY = 0;

    int rank_esq = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    int rank_dir = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

    srand(42 + rank); 
    
    std::vector<std::vector<Cell>> grid_local(local_W, std::vector<Cell>(local_H));
    for (int i = 0; i < local_W; ++i) {
        for (int j = 0; j < local_H; ++j) {
            int gx = offsetX + i;
            int gy = offsetY + j;
            
            int bloco_x = gx / 40; 
            int bloco_y = gy / 40;
            grid_local[i][j].tipo = (bloco_x + bloco_y * 3) % 5; 
            
            grid_local[i][j].recurso = (grid_local[i][j].tipo != 4) ? (rand() % 40) : 0; 
            
            grid_local[i][j].acessivel = (grid_local[i][j].tipo != 4); 
            grid_local[i][j].consumo_acumulado = 0;
        }
    }

    std::vector<Agent> agentes_locais;
    
    int agentes_por_processo = N_AGENTS / size;
    int resto_agentes = N_AGENTS % size;

    int meus_agentes = agentes_por_processo + (rank == 0 ? resto_agentes : 0);

    int id_start = rank * agentes_por_processo + (rank > 0 ? resto_agentes : 0);

    for (int i = 0; i < meus_agentes; ++i) {
        Agent a;
        a.id = id_start + i;
        
        int lx, ly;
        do {
            lx = rand() % local_W;
            ly = rand() % local_H;
        } while (grid_local[lx][ly].acessivel == false);
        
        a.x = offsetX + lx;
        a.y = offsetY + ly;
        
        a.energia = 100;
        
        agentes_locais.push_back(a);
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    for(int p = 0; p < size; ++p) {
        if(rank == p) {
            std::cout << "[INIT] Processo " << rank << " iniciou com subgrid " 
                      << local_W << "x" << local_H << " e " 
                      << agentes_locais.size() << " agentes." << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int estacao = 0;

    for (int t = 0; t < T; ++t) {
        
        if (t % S == 0) {
            if (rank == 0) {
                estacao = 1 - estacao;
            }
            MPI_Bcast(&estacao, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }

        std::vector<Cell> halo_recebido_esq(local_H);
        std::vector<Cell> halo_recebido_dir(local_H);
        MPI_Request reqs_halo[4];

        MPI_Isend(grid_local[0].data(), local_H, MPI_CELL, rank_esq, 0, MPI_COMM_WORLD, &reqs_halo[0]);
        MPI_Irecv(halo_recebido_dir.data(), local_H, MPI_CELL, rank_dir, 0, MPI_COMM_WORLD, &reqs_halo[1]);

        MPI_Isend(grid_local[local_W - 1].data(), local_H, MPI_CELL, rank_dir, 1, MPI_COMM_WORLD, &reqs_halo[2]);
        MPI_Irecv(halo_recebido_esq.data(), local_H, MPI_CELL, rank_esq, 1, MPI_COMM_WORLD, &reqs_halo[3]);

        MPI_Waitall(4, reqs_halo, MPI_STATUSES_IGNORE);

        std::vector<Agent> buffer_envio_esquerda, buffer_envio_direita;
        std::vector<Agent> nova_lista_local;

        #pragma omp parallel 
        {
            std::vector<Agent> lista_local_thread; 
            std::vector<Agent> envio_esq_thread, envio_dir_thread; 

            #pragma omp for 
            for (size_t k = 0; k < agentes_locais.size(); ++k) { 
                Agent a = agentes_locais[k];
                int lx = a.x - offsetX;
                int r = grid_local[lx][a.y].recurso; 

                executar_carga(r); 

                Destino dest = decidir(a, offsetX, local_W, local_H, 
                                       grid_local, halo_recebido_esq, halo_recebido_dir, 
                                       rank, rank_esq, rank_dir);

                a.x = dest.novo_x;
                a.y = dest.novo_y;

                if (dest.eh_local) {
                    int novo_lx = a.x - offsetX;
                    
                    #pragma omp atomic
                    grid_local[novo_lx][a.y].consumo_acumulado += 1; 
                    
                    lista_local_thread.push_back(a); 
                } else {
                    if (dest.processo_destino == rank_esq) {
                        envio_esq_thread.push_back(a);
                    } else if (dest.processo_destino == rank_dir) {
                        envio_dir_thread.push_back(a);
                    }
                }
            }

            #pragma omp critical
            {
                nova_lista_local.insert(nova_lista_local.end(), lista_local_thread.begin(), lista_local_thread.end());
                buffer_envio_esquerda.insert(buffer_envio_esquerda.end(), envio_esq_thread.begin(), envio_esq_thread.end());
                buffer_envio_direita.insert(buffer_envio_direita.end(), envio_dir_thread.begin(), envio_dir_thread.end());
            }
        }
        
        agentes_locais = nova_lista_local;

        int qtd_envio_esq = buffer_envio_esquerda.size();
        int qtd_envio_dir = buffer_envio_direita.size();
        int qtd_recv_esq = 0, qtd_recv_dir = 0;

        MPI_Request reqs_qtd[4];
        MPI_Isend(&qtd_envio_esq, 1, MPI_INT, rank_esq, 2, MPI_COMM_WORLD, &reqs_qtd[0]);
        MPI_Irecv(&qtd_recv_dir, 1, MPI_INT, rank_dir, 2, MPI_COMM_WORLD, &reqs_qtd[1]);
        MPI_Isend(&qtd_envio_dir, 1, MPI_INT, rank_dir, 3, MPI_COMM_WORLD, &reqs_qtd[2]);
        MPI_Irecv(&qtd_recv_esq, 1, MPI_INT, rank_esq, 3, MPI_COMM_WORLD, &reqs_qtd[3]);
        MPI_Waitall(4, reqs_qtd, MPI_STATUSES_IGNORE);

        std::vector<Agent> buffer_recv_esq(qtd_recv_esq);
        std::vector<Agent> buffer_recv_dir(qtd_recv_dir);
        MPI_Request reqs_agentes[4];

        MPI_Isend(buffer_envio_esquerda.data(), qtd_envio_esq, MPI_AGENT, rank_esq, 4, MPI_COMM_WORLD, &reqs_agentes[0]);
        MPI_Irecv(buffer_recv_dir.data(), qtd_recv_dir, MPI_AGENT, rank_dir, 4, MPI_COMM_WORLD, &reqs_agentes[1]);
        MPI_Isend(buffer_envio_direita.data(), qtd_envio_dir, MPI_AGENT, rank_dir, 5, MPI_COMM_WORLD, &reqs_agentes[2]);
        MPI_Irecv(buffer_recv_esq.data(), qtd_recv_esq, MPI_AGENT, rank_esq, 5, MPI_COMM_WORLD, &reqs_agentes[3]);
        MPI_Waitall(4, reqs_agentes, MPI_STATUSES_IGNORE);

        agentes_locais.insert(agentes_locais.end(), buffer_recv_esq.begin(), buffer_recv_esq.end());
        agentes_locais.insert(agentes_locais.end(), buffer_recv_dir.begin(), buffer_recv_dir.end());

        std::string arquivo_agentes = "agentes_rank_" + std::to_string(rank) + ".csv";
        std::ofstream file_agentes(arquivo_agentes, std::ios_base::app);
        if (t == 0) file_agentes << "ciclo,id,x,y\n"; 
        for (size_t k = 0; k < agentes_locais.size(); ++k) {
            file_agentes << t << "," << agentes_locais[k].id << "," 
                         << agentes_locais[k].x << "," << agentes_locais[k].y << "\n";
        }
        file_agentes.close();

        std::string arquivo_grid = "grid_rank_" + std::to_string(rank) + ".csv";
        std::ofstream file_grid(arquivo_grid, std::ios_base::app);
        if (t == 0) file_grid << "ciclo,x,y,recurso\n"; 
        for (int i = 0; i < local_W; ++i) {
            for (int j = 0; j < local_H; ++j) {
                file_grid << t << "," << (offsetX + i) << "," << (offsetY + j) << "," 
                          << grid_local[i][j].recurso << "\n";
            }
        }
        file_grid.close();

        #pragma omp parallel for collapse(2) 
        for (int i = 0; i < local_W; ++i) { 
            for (int j = 0; j < local_H; ++j) {
                
                int tipo = grid_local[i][j].tipo;
                int taxa_regen = 0;

                if (tipo == 0) {
                    taxa_regen = 1;
                } else if (tipo == 1) {
                    taxa_regen = (estacao == 1) ? 1 : 0; 
                } else if (tipo == 2) {
                    taxa_regen = 1; 
                } else if (tipo == 3) {
                    taxa_regen = (estacao == 0) ? 1 : 0; 
                } else if (tipo == 4) {
                    taxa_regen = 0;
                }

                grid_local[i][j].recurso += taxa_regen;  
                
                grid_local[i][j].recurso -= (grid_local[i][j].consumo_acumulado * 50); 
                
                if(grid_local[i][j].recurso < 0) grid_local[i][j].recurso = 0;
                if(grid_local[i][j].recurso > 100) grid_local[i][j].recurso = 100;
                
                grid_local[i][j].consumo_acumulado = 0; 
            }
        }

        int total_agentes_locais = agentes_locais.size();
        int total_agentes_globais = 0;
        MPI_Allreduce(&total_agentes_locais, &total_agentes_globais, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        long long recurso_local_total = 0;
        
        #pragma omp parallel for reduction(+:recurso_local_total) collapse(2)
        for (int i = 0; i < local_W; ++i) {
            for (int j = 0; j < local_H; ++j) {
                recurso_local_total += grid_local[i][j].recurso;
            }
        }

        long long recurso_global_total = 0;
        MPI_Reduce(&recurso_local_total, &recurso_global_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            if (t % 50 == 0 || t == T - 1) { 
                std::string nome_estacao = (estacao == 0) ? "Seca" : "Cheia";
                std::cout << "[Ciclo " << t << "] Estação: " << nome_estacao 
                          << " | Agentes Globais: " << total_agentes_globais 
                          << " | Recurso Total no Território: " << recurso_global_total
                          << "\n" << std::endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0) {
        std::cout << "---" << std::endl;
        std::cout << "Simulacao finalizada com sucesso em " << MPI_Wtime() - start_time << " segundos." << std::endl;
    }

    MPI_Type_free(&MPI_AGENT);
    MPI_Type_free(&MPI_CELL);

    MPI_Finalize();
    return 0;
}