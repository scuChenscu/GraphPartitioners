#pragma once

// 节点结构
template <typename K, typename V>
struct Node {
    K key;
    V value;
    Node* next;

    Node(K k, V v) : key(k), value(v), next(nullptr) {}
};

// 链表类
template <typename K, typename V>
class LinkedList {
private:
    Node<K, V>* head;

public:
    LinkedList() : head(nullptr) {}

    // 判空函数
    bool isEmpty() const {
        return head == nullptr;
    }

    // 插入键值对
    void insert(K key, V value) {
        Node<K, V>* newNode = new Node<K, V>(key, value);
        newNode->next = head;
        head = newNode;
    }

    // 移除指定 key 的元素
    std::pair<K, V> removeByKey(K key) {
        if (isEmpty()) {
            std::cerr << "Linked List is empty." << std::endl;
            return {K(), V()};  // 返回默认键值对
        }

        Node<K, V>* current = head;
        Node<K, V>* prev = nullptr;

        while (current != nullptr) {
            if (current->key == key) {
                if (prev == nullptr) {
                    // 如果是头结点
                    head = current->next;
                } else {
                    // 如果是中间或尾部节点
                    prev->next = current->next;
                }

                std::pair<K, V> removedPair = {current->key, current->value};
                delete current;
                return removedPair;
            }

            prev = current;
            current = current->next;
        }

        std::cerr << "Element with key '" << key << "' not found." << std::endl;
        return {K(), V()};  // 返回默认键值对
    }

    // 移除首个元素
    std::pair<K, V> removeFirst() {
        if (isEmpty()) {
            std::cerr << "Linked List is empty." << std::endl;
            return {K(), V()};  // 返回默认键值对
        }

        Node<K, V>* temp = head;
        head = head->next;
        std::pair<K, V> removedPair = {temp->key, temp->value};
        delete temp;
        return removedPair;
    }

    // 修改指定 key 的值
    void modifyValue(K key, V newValue) {
        if (isEmpty()) {
            std::cerr << "Linked List is empty." << std::endl;
            return;
        }

        Node<K, V>* current = head;

        while (current != nullptr) {
            if (current->key == key) {
                current->value = newValue;
                return;
            }

            current = current->next;
        }

        std::cerr << "Element with key '" << key << "' not found." << std::endl;
    }

    // 根据 key 查询元素的 value
    V getValueByKey(K key) const {
        if (isEmpty()) {
            std::cerr << "Linked List is empty." << std::endl;
            return V();  // 返回默认值
        }

        Node<K, V>* current = head;

        while (current != nullptr) {
            if (current->key == key) {
                return current->value;
            }

            current = current->next;
        }

        std::cerr << "Element with key '" << key << "' not found." << std::endl;
        return V();  // 返回默认值
    }

    // 打印链表内容
    void printList() {
        if (isEmpty()) {
            std::cout << "Linked List is empty." << std::endl;
            return;
        }

        Node<K, V>* current = head;

        while (current != nullptr) {
            std::cout << "(" << current->key << ", " << current->value << ") ";
            current = current->next;
        }

        std::cout << std::endl;
    }
};

//int main() {
//    LinkedList<int, std::string> linkedList;
//
//    // 判空
//    std::cout << "Is Linked List empty? " << (linkedList.isEmpty() ? "Yes" : "No") << std::endl;
//
//    // 插入键值对
//    linkedList.insert(1, "One");
//    linkedList.insert(2, "Two");
//    linkedList.insert(3, "Three");
//
//    // 打印链表内容
//    std::cout << "Original Linked List: ";
//    linkedList.printList();
//
//    // 查询指定 key 的元素值
//    int keyToQuery = 2;
//    std::string value = linkedList.getValueByKey(keyToQuery);
//    std::cout << "Value for key " << keyToQuery << ": " << value << std::endl;
//
//    // 判空
//    std::cout << "Is Linked List empty? " << (linkedList.isEmpty() ? "Yes" : "No") << std::endl;
//
//    return 0;
//}
