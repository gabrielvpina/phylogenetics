import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction # CORRIGIDO: Agora usa gc_fraction

def fasta_para_tsv(arquivo_fasta):
    """
    Lê um arquivo FASTA e imprime o nome, tamanho e o conteúdo GC (%) 
    de cada sequência no formato TSV (separado por tabulação) no terminal.
    """
    if not arquivo_fasta:
        print("Uso: python fasta_info.py <arquivo.fasta>", file=sys.stderr)
        sys.exit(1)

    try:
        # Imprime o novo cabeçalho da tabela, incluindo a coluna GC
        print("Nome_da_Sequencia\tTamanho\tConteudo_GC_%")
        
        # O SeqIO lida com múltiplas sequências e quebras de linha no FASTA
        for registro in SeqIO.parse(arquivo_fasta, "fasta"):
            
            sequencia = registro.seq
            tamanho = len(sequencia)
            
            # 1. Calcula o conteúdo GC usando a função CORRIGIDA
            try:
                # gc_fraction retorna um valor entre 0 e 1. Multiplicamos por 100.
                gc_valor = gc_fraction(sequencia) * 100
                gc_porcentagem = round(gc_valor, 2)
            except ValueError:
                # Trata sequências que não são DNA/RNA (e.g., só contém N's)
                gc_porcentagem = "N/A" 
                
            # 2. Imprime os dados no formato TSV (separados por \t)
            print(f"{registro.id}\t{tamanho}\t{gc_porcentagem}")
            
    except FileNotFoundError:
        print(f"Erro: Arquivo '{arquivo_fasta}' não encontrado.", file=sys.stderr)
        sys.exit(1)
    except ImportError:
        # Verifica se o Biopython está instalado, necessário para SeqUtils
        print("Erro: A biblioteca 'Biopython' não está instalada. Instale com 'pip install biopython'.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Ocorreu um erro: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # O script espera o caminho do arquivo FASTA como o primeiro argumento
    arquivo_fasta = sys.argv[1] if len(sys.argv) > 1 else None
    fasta_para_tsv(arquivo_fasta)
